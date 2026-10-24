%% Part C - relay/destination reception and measured metrics (base MATLAB)
% Put this file beside partA_data.mat and partB_data.mat, then run it.
% The source/relay framing and channel observations are from Parts A and B.
clear; clc; close all;
outputFolder=fileparts(mfilename('fullpath'));
fileA=fullfile(outputFolder,'partA_data.mat');
fileB=fullfile(outputFolder,'partB_data.mat');
assert(isfile(fileA) && isfile(fileB), ...
    'Run Parts A and B first; keep both .mat files beside this script.');
A=load(fileA,'tx'); B=load(fileB);
tx=A.tx;
positions=B.cfg.relayPositionsM;
EbN0dB=B.cfg.EbN0dB;
nPos=numel(positions);
nSNR=numel(EbN0dB);
nFrames=tx.cfg.numberOfFrames;
nSym=tx.symbolsPerFrame;
nPayload=numel(tx.payloadRows);
assert(nFrames==size(B.rxSD,2) && nSym==size(B.rxSD,1));
assert(isequal(size(B.rxSR),[nSym nFrames nSNR nPos]));
assert(all(B.N0>0) && numel(B.N0)==nSNR);
assert(tx.crcBits==8 && nPayload==tx.cfg.payloadBitsPerFrame);
for f=1:nFrames
    assert(crc_ok_c(tx.sourceFrames(:,f)), ...
        'Part A CRC does not match the receiver CRC implementation.');
end
assert(max(abs(tx.sourceEnergy/nSym-tx.sourcePower))<1e-12);
assert(max(abs(tx.relayEnergy/nSym-tx.relayPower))<1e-12);

% Outage is explicitly a link-symbol-SNR event. The same 5 dB threshold
% applies to direct reception and to the selected/combined cooperative link.
outageThresholdDB=5;
outageThreshold=10^(outageThresholdDB/10);
BER_direct=zeros(1,nSNR);
PER_direct=zeros(1,nSNR);
OUT_direct=zeros(1,nSNR);
TH_direct=zeros(1,nSNR);       % useful bits per symbol-period time slot
BER_coop=zeros(nPos,nSNR);
PER_coop=zeros(nPos,nSNR);
OUT_coop=zeros(nPos,nSNR);
TH_coop=zeros(nPos,nSNR);
relayForwardRate=zeros(nPos,nSNR);
CRCfail_direct=zeros(1,nSNR);
CRCfail_coop=zeros(nPos,nSNR);
undetected_direct=zeros(1,nSNR);
undetected_coop=zeros(nPos,nSNR);
relayCRCfail=zeros(nPos,nSNR);

%% Direct one-slot receiver: full frame energy, same S-D fading/noise as
% the source observation in the cooperative two-slot comparison.
directAmplitude=sqrt(B.gainSD*B.EsReference);
for k=1:nSNR
    bitErrors=0; packetErrors=0; outages=0; successes=0;
    crcFailures=0; undetected=0;
    for f=1:nFrames
        aD=directAmplitude*B.hSD(f,k);
        estimate=B.rxDirect(:,f,k)/aD;
        decoded=qpsk_demap_c(estimate);
        truth=tx.sourceFrames(:,f);
        badPacket=any(decoded~=truth);
        crcOK=crc_ok_c(decoded);
        bitErrors=bitErrors+nnz(decoded(tx.payloadRows)~=tx.payloadBits(:,f));
        packetErrors=packetErrors+badPacket;
        crcFailures=crcFailures+~crcOK;
        undetected=undetected+(crcOK && badPacket);
        successes=successes+(crcOK && ~badPacket);
        outages=outages+(abs(aD)^2/B.N0(k)<outageThreshold);
    end
    BER_direct(k)=bitErrors/(nFrames*nPayload);
    PER_direct(k)=packetErrors/nFrames;
    OUT_direct(k)=outages/nFrames;
    TH_direct(k)=successes*nPayload/(nFrames*nSym);
    CRCfail_direct(k)=crcFailures/nFrames;
    undetected_direct(k)=undetected;
end

%% Decode-and-forward relay; MRC of two independently faded observations.
% In phase 1, both the relay and destination listen to the same source
% frame. In phase 2, the relay transmits only after its CRC passes.
for p=1:nPos
    for k=1:nSNR
        bitErrors=0; packetErrors=0; outages=0; successes=0;
        forwards=0; relayFailures=0; crcFailures=0; undetected=0;
        for f=1:nFrames
            truth=tx.sourceFrames(:,f);
            aSR=sqrt(B.gainSR(p)*tx.sourcePower)*B.hSR(f,k,p);
            decodedRelay=qpsk_demap_c(B.rxSR(:,f,k,p)/aSR);
            relayOK=crc_ok_c(decodedRelay);
            forwards=forwards+relayOK;
            relayFailures=relayFailures+~relayOK;

            aSD=sqrt(B.gainSD*tx.sourcePower)*B.hSD(f,k);
            ySD=B.rxSD(:,f,k);
            if relayOK
                relayTx=sqrt(tx.relayPower)*qpsk_map_c(decodedRelay);
                aRD=sqrt(B.gainRD(p)*tx.relayPower)*B.hRD(f,k,p);
                yRD=sqrt(B.gainRD(p))*B.hRD(f,k,p)*relayTx + ...
                    B.noiseRD(:,f,k,p);
                % Equal-noise-variance MRC; use both observations of the
                % same QPSK frame. A rare undetected relay CRC error is
                % tracked rather than assumed impossible.
                combined=(conj(aSD)*ySD+conj(aRD)*yRD)/ ...
                    (abs(aSD)^2+abs(aRD)^2);
                effectiveSNR=(abs(aSD)^2+abs(aRD)^2)/B.N0(k);
            else
                % Relay is silent: do not combine R-D receiver noise.
                combined=ySD/aSD;
                effectiveSNR=abs(aSD)^2/B.N0(k);
            end
            decoded=qpsk_demap_c(combined);
            badPacket=any(decoded~=truth);
            crcOK=crc_ok_c(decoded);
            bitErrors=bitErrors+nnz(decoded(tx.payloadRows)~=tx.payloadBits(:,f));
            packetErrors=packetErrors+badPacket;
            crcFailures=crcFailures+~crcOK;
            undetected=undetected+(crcOK && badPacket);
            successes=successes+(crcOK && ~badPacket);
            outages=outages+(effectiveSNR<outageThreshold);
        end
        BER_coop(p,k)=bitErrors/(nFrames*nPayload);
        PER_coop(p,k)=packetErrors/nFrames;
        OUT_coop(p,k)=outages/nFrames;
        % Both time slots are reserved per cooperative attempt, including
        % a silent relay. Each time slot contains nSym symbol periods.
        TH_coop(p,k)=successes*nPayload/(nFrames*2*nSym);
        relayForwardRate(p,k)=forwards/nFrames;
        relayCRCfail(p,k)=relayFailures/nFrames;
        CRCfail_coop(p,k)=crcFailures/nFrames;
        undetected_coop(p,k)=undetected;
    end
end

assert(all(BER_direct>=0 & BER_direct<=1));
assert(all(BER_coop(:)>=0 & BER_coop(:)<=1));
assert(all(TH_direct<=nPayload/nSym+1e-12));
assert(all(TH_coop(:)<=nPayload/(2*nSym)+1e-12));
assert(all(abs(relayForwardRate(:)+relayCRCfail(:)-1)<1e-12));

%% Export actual numerical results (never substitute plot floors in CSV).
PositionM=[]; NominalEbN0dB=[]; Mode={}; BER=[]; PER=[];
Outage=[]; UsefulBitsPerSymbolPeriod=[]; RelayForwardRate=[];
for k=1:nSNR
    PositionM(end+1,1)=NaN; %#ok<SAGROW> direct has no relay position
    NominalEbN0dB(end+1,1)=EbN0dB(k); %#ok<SAGROW>
    Mode{end+1,1}='Direct'; %#ok<SAGROW>
    BER(end+1,1)=BER_direct(k); %#ok<SAGROW>
    PER(end+1,1)=PER_direct(k); %#ok<SAGROW>
    Outage(end+1,1)=OUT_direct(k); %#ok<SAGROW>
    UsefulBitsPerSymbolPeriod(end+1,1)=TH_direct(k); %#ok<SAGROW>
    RelayForwardRate(end+1,1)=NaN; %#ok<SAGROW>
    for p=1:nPos
        PositionM(end+1,1)=positions(p); %#ok<SAGROW>
        NominalEbN0dB(end+1,1)=EbN0dB(k); %#ok<SAGROW>
        Mode{end+1,1}='Cooperative'; %#ok<SAGROW>
        BER(end+1,1)=BER_coop(p,k); %#ok<SAGROW>
        PER(end+1,1)=PER_coop(p,k); %#ok<SAGROW>
        Outage(end+1,1)=OUT_coop(p,k); %#ok<SAGROW>
        UsefulBitsPerSymbolPeriod(end+1,1)=TH_coop(p,k); %#ok<SAGROW>
        RelayForwardRate(end+1,1)=relayForwardRate(p,k); %#ok<SAGROW>
    end
end
results=table(Mode,PositionM,NominalEbN0dB,BER,PER,Outage, ...
    UsefulBitsPerSymbolPeriod,RelayForwardRate);
writetable(results,fullfile(outputFolder,'partC_results.csv'));
save(fullfile(outputFolder,'partC_data.mat'),'EbN0dB','positions', ...
    'outageThresholdDB','BER_direct','PER_direct','OUT_direct', ...
    'TH_direct','BER_coop','PER_coop','OUT_coop','TH_coop', ...
    'relayForwardRate','CRCfail_direct','CRCfail_coop', ...
    'relayCRCfail','undetected_direct','undetected_coop','results');

%% Evidence plots. A floor is used for log-axis display only; the CSV
% retains exact zeros and a zero means no errors observed in this sample.
berFloor=1/(nFrames*nPayload);
perFloor=1/nFrames;
labels=[{'Direct'},arrayfun(@(x)sprintf('Cooperative %g m',x), ...
    positions,'UniformOutput',false)];
colors=lines(nPos);
figure('Color','w');
semilogy(EbN0dB,max(BER_direct,berFloor),'k-o','LineWidth',1.5);
hold on;
for p=1:nPos
    semilogy(EbN0dB,max(BER_coop(p,:),berFloor),'-s', ...
        'Color',colors(p,:),'LineWidth',1.4);
end
grid on; xlabel('Nominal direct-reference E_b/N_0 (dB)');
ylabel('Payload BER'); title('Direct versus CRC-gated cooperative BER');
legend(labels,'Location','southwest');
exportgraphics(gcf,fullfile(outputFolder,'partC_BER.png'),'Resolution',200);

figure('Color','w');
semilogy(EbN0dB,max(PER_direct,perFloor),'k-o','LineWidth',1.5);
hold on;
for p=1:nPos
    semilogy(EbN0dB,max(PER_coop(p,:),perFloor),'-s', ...
        'Color',colors(p,:),'LineWidth',1.4);
end
grid on; xlabel('Nominal direct-reference E_b/N_0 (dB)');
ylabel('True frame error rate (PER)');
title('Direct versus CRC-gated cooperative packet errors');
legend(labels,'Location','southwest');
exportgraphics(gcf,fullfile(outputFolder,'partC_PER.png'),'Resolution',200);

figure('Color','w');
plot(EbN0dB,OUT_direct,'k-o','LineWidth',1.5); hold on;
for p=1:nPos
    plot(EbN0dB,OUT_coop(p,:),'-s', ...
        'Color',colors(p,:),'LineWidth',1.4);
end
grid on; ylim([0 1]); xlabel('Nominal direct-reference E_b/N_0 (dB)');
ylabel('Outage probability');
title(sprintf('Outage: instantaneous symbol SNR < %g dB',outageThresholdDB));
legend(labels,'Location','northeast');
exportgraphics(gcf,fullfile(outputFolder,'partC_outage.png'),'Resolution',200);

figure('Color','w');
plot(EbN0dB,TH_direct,'k-o','LineWidth',1.5); hold on;
for p=1:nPos
    plot(EbN0dB,TH_coop(p,:),'-s', ...
        'Color',colors(p,:),'LineWidth',1.4);
end
grid on; xlabel('Nominal direct-reference E_b/N_0 (dB)');
ylabel('Delivered payload bits per symbol period');
title('Effective throughput with CRC and second-slot overhead');
legend(labels,'Location','southeast');
exportgraphics(gcf,fullfile(outputFolder,'partC_throughput.png'),'Resolution',200);

figure('Color','w'); hold on;
for p=1:nPos
    plot(EbN0dB,relayForwardRate(p,:),'-o', ...
        'Color',colors(p,:),'LineWidth',1.4);
end
grid on; ylim([0 1]); xlabel('Nominal direct-reference E_b/N_0 (dB)');
ylabel('Fraction of frames forwarded');
title('Relay forwarding after S-R CRC validation');
legend(arrayfun(@(x)sprintf('%g m',x),positions, ...
    'UniformOutput',false),'Location','southeast');
exportgraphics(gcf,fullfile(outputFolder,'partC_relay_forwarding.png'), ...
    'Resolution',200);

fprintf('PART C COMPLETE: %d frames x %d Eb/N0 values x %d positions.\n', ...
    nFrames,nSNR,nPos);
fprintf('Outage threshold: %g dB instantaneous symbol SNR.\n', ...
    outageThresholdDB);
fprintf('At %g dB nominal reference:\n',EbN0dB(end));
fprintf('  Direct: BER %.6g, PER %.4f, outage %.4f, throughput %.4f\n', ...
    BER_direct(end),PER_direct(end),OUT_direct(end),TH_direct(end));
for p=1:nPos
    fprintf('  Relay %g m: BER %.6g, PER %.4f, outage %.4f, throughput %.4f, forwarding %.3f\n', ...
        positions(p),BER_coop(p,end),PER_coop(p,end),OUT_coop(p,end), ...
        TH_coop(p,end),relayForwardRate(p,end));
end
fprintf('Numerical results: partC_results.csv; plots: partC_*.png\n');

function bits=qpsk_demap_c(symbols)
bits=zeros(2*numel(symbols),1);
bits(1:2:end)=real(symbols(:))<0;
bits(2:2:end)=imag(symbols(:))<0;
end

function symbols=qpsk_map_c(bits)
bits=bits(:);
assert(mod(numel(bits),2)==0);
symbols=((1-2*bits(1:2:end))+1i*(1-2*bits(2:2:end)))/sqrt(2);
end

function valid=crc_ok_c(frame)
frame=frame(:);
valid=isequal(frame(end-7:end),crc8_bits_c(frame(1:end-8)));
end

function remainder=crc8_bits_c(bits)
% Same CRC-8: g(x)=x^8+x^2+x+1, init 0, no reflection/final XOR.
bits=double(bits(:));
work=[bits;zeros(8,1)];
generator=[1 0 0 0 0 0 1 1 1]';
for i=1:numel(bits)
    if work(i)==1
        work(i:i+8)=mod(work(i:i+8)+generator,2);
    end
end
remainder=work(end-7:end);
end

%% Part D - engineering comparison and relay placement (base MATLAB)
% Run Parts A, B and C first, then place this script beside their .mat files.
% It reuses the exact same fading/noise samples for all power allocations.
clear; clc; close all;
outputFolder=fileparts(mfilename('fullpath'));
fileA=fullfile(outputFolder,'partA_data.mat');
fileB=fullfile(outputFolder,'partB_data.mat');
fileC=fullfile(outputFolder,'partC_data.mat');
assert(isfile(fileA) && isfile(fileB) && isfile(fileC), ...
    'Run Parts A, B and C first; keep their data files beside Part D.');
A=load(fileA,'tx'); B=load(fileB); C=load(fileC);
tx=A.tx;
positions=B.cfg.relayPositionsM;
EbN0dB=B.cfg.EbN0dB;
nPos=numel(positions); nSNR=numel(EbN0dB);
nFrames=tx.cfg.numberOfFrames; nSym=tx.symbolsPerFrame;
nPayload=numel(tx.payloadRows);
assert(isequal(C.positions,positions) && isequal(C.EbN0dB,EbN0dB));
assert(nFrames==size(B.rxSR,2) && nSym==size(B.rxSR,1));
assert(tx.crcBits==8 && nPayload==tx.cfg.payloadBitsPerFrame);

%% Stated engineering target - editable before running.
targetEbN0dB=12;              % nominal full-budget direct reference
targetBER=0.003;              % payload bit error rate <= 0.3%
targetOutage=0.05;            % instantaneous symbol-SNR outage <= 5%
outageThresholdDB=C.outageThresholdDB; % 5 dB from Part C
threshold=10^(outageThresholdDB/10);
targetIndex=find(EbN0dB==targetEbN0dB,1);
assert(~isempty(targetIndex),'Target Eb/N0 must be in the Part B sweep.');

% This is an energy-per-attempt cap, not a guarantee that all 500 units
% are consumed. If the relay CRC fails, phase 2 is silent and its energy
% remains unused; a second time slot is still reserved in the comparison.
sourceFraction=[0.25 0.50 0.75];
nSplit=numel(sourceFraction);
totalEnergy=tx.cfg.totalEnergyPerFrame;
Psource=sourceFraction*totalEnergy/nSym;
Prelay=(1-sourceFraction)*totalEnergy/nSym;
assert(all(abs((Psource+Prelay)*nSym-totalEnergy)<1e-9));

BER=zeros(nSplit,nPos,nSNR);
PER=zeros(nSplit,nPos,nSNR);
OUT=zeros(nSplit,nPos,nSNR);
TH=zeros(nSplit,nPos,nSNR);
forwardRate=zeros(nSplit,nPos,nSNR);
undetectedRelay=zeros(nSplit,nPos,nSNR);

%% Rerun reception under changed powers on the SAME noise/fading draws.
% Part B saved observations at the Part A default source power. Subtract
% their exact noiseless component to recover each noise realization, then
% add the source signal at the new power. No extra random numbers are used.
for a=1:nSplit
    for p=1:nPos
        for k=1:nSNR
            bitErrors=0; packetErrors=0; outages=0;
            delivered=0; forwarded=0; relayUndetected=0;
            for f=1:nFrames
                x=tx.sourceSymbols(:,f);
                truth=tx.sourceFrames(:,f);
                hSD=B.hSD(f,k); hSR=B.hSR(f,k,p);
                hRD=B.hRD(f,k,p);
                aSDold=sqrt(B.gainSD*tx.sourcePower)*hSD;
                aSRold=sqrt(B.gainSR(p)*tx.sourcePower)*hSR;
                nSD=B.rxSD(:,f,k)-aSDold*x;
                nSR=B.rxSR(:,f,k,p)-aSRold*x;

                aSD=sqrt(B.gainSD*Psource(a))*hSD;
                aSR=sqrt(B.gainSR(p)*Psource(a))*hSR;
                ySD=aSD*x+nSD;
                ySR=aSR*x+nSR;
                bitsRelay=qpsk_demap_d(ySR/aSR);
                relayOK=crc_ok_d(bitsRelay);
                forwarded=forwarded+relayOK;
                relayUndetected=relayUndetected+ ...
                    (relayOK && any(bitsRelay~=truth));

                if relayOK
                    aRD=sqrt(B.gainRD(p)*Prelay(a))*hRD;
                    yRD=aRD*qpsk_map_d(bitsRelay)+B.noiseRD(:,f,k,p);
                    equalized=(conj(aSD)*ySD+conj(aRD)*yRD)/ ...
                        (abs(aSD)^2+abs(aRD)^2);
                    effectiveSNR=(abs(aSD)^2+abs(aRD)^2)/B.N0(k);
                else
                    % Silent relay: use the direct source observation only.
                    equalized=ySD/aSD;
                    effectiveSNR=abs(aSD)^2/B.N0(k);
                end

                decoded=qpsk_demap_d(equalized);
                correct=isequal(decoded,truth);
                destinationCRC=crc_ok_d(decoded);
                bitErrors=bitErrors+nnz(decoded(tx.payloadRows)~= ...
                    tx.payloadBits(:,f));
                packetErrors=packetErrors+~correct;
                delivered=delivered+(correct && destinationCRC);
                outages=outages+(effectiveSNR<threshold);
            end
            BER(a,p,k)=bitErrors/(nFrames*nPayload);
            PER(a,p,k)=packetErrors/nFrames;
            OUT(a,p,k)=outages/nFrames;
            TH(a,p,k)=delivered*nPayload/(nFrames*2*nSym);
            forwardRate(a,p,k)=forwarded/nFrames;
            undetectedRelay(a,p,k)=relayUndetected;
        end
    end
end

%% Cross-check the 50/50 rerun against Part C before drawing conclusions.
defaultIndex=find(abs(sourceFraction-tx.cfg.sourceEnergyFraction)<1e-12,1);
assert(~isempty(defaultIndex),'Include the Part A split in Part D.');
maxBERDiff=max(abs(reshape(BER(defaultIndex,:,:),nPos,nSNR)-C.BER_coop),[],'all');
maxPERDiff=max(abs(reshape(PER(defaultIndex,:,:),nPos,nSNR)-C.PER_coop),[],'all');
maxOUTDiff=max(abs(reshape(OUT(defaultIndex,:,:),nPos,nSNR)-C.OUT_coop),[],'all');
maxTHDiff=max(abs(reshape(TH(defaultIndex,:,:),nPos,nSNR)-C.TH_coop),[],'all');
assert(maxBERDiff<1/(nFrames*nPayload)+1e-12 && ...
    maxPERDiff<1/nFrames+1e-12 && maxOUTDiff<1/nFrames+1e-12 && ...
    maxTHDiff<nPayload/(nFrames*2*nSym)+1e-12, ...
    'The 50/50 rerun differs from Part C; check input files and power units.');

%% Full table and per-position improvement regions.
Alpha=[]; PositionM=[]; NominalEbN0dB=[]; CoopBER=[]; DirectBER=[];
CoopPER=[]; DirectPER=[]; CoopOutage=[]; DirectOutage=[];
CoopThroughput=[]; DirectThroughput=[]; ForwardRate=[];
for a=1:nSplit
    for p=1:nPos
        for k=1:nSNR
            Alpha(end+1,1)=sourceFraction(a); %#ok<SAGROW>
            PositionM(end+1,1)=positions(p); %#ok<SAGROW>
            NominalEbN0dB(end+1,1)=EbN0dB(k); %#ok<SAGROW>
            CoopBER(end+1,1)=BER(a,p,k); %#ok<SAGROW>
            DirectBER(end+1,1)=C.BER_direct(k); %#ok<SAGROW>
            CoopPER(end+1,1)=PER(a,p,k); %#ok<SAGROW>
            DirectPER(end+1,1)=C.PER_direct(k); %#ok<SAGROW>
            CoopOutage(end+1,1)=OUT(a,p,k); %#ok<SAGROW>
            DirectOutage(end+1,1)=C.OUT_direct(k); %#ok<SAGROW>
            CoopThroughput(end+1,1)=TH(a,p,k); %#ok<SAGROW>
            DirectThroughput(end+1,1)=C.TH_direct(k); %#ok<SAGROW>
            ForwardRate(end+1,1)=forwardRate(a,p,k); %#ok<SAGROW>
        end
    end
end
ImprovesBER=CoopBER<DirectBER;
ImprovesPER=CoopPER<DirectPER;
ImprovesOutage=CoopOutage<DirectOutage;
LosesThroughput=CoopThroughput<DirectThroughput;
allResults=table(Alpha,PositionM,NominalEbN0dB,CoopBER,DirectBER, ...
    CoopPER,DirectPER,CoopOutage,DirectOutage,CoopThroughput, ...
    DirectThroughput,ForwardRate,ImprovesBER,ImprovesPER, ...
    ImprovesOutage,LosesThroughput);
writetable(allResults,fullfile(outputFolder,'partD_all_results.csv'));

%% Select a design ONLY from those meeting both stated targets at 12 dB.
targetMask=NominalEbN0dB==targetEbN0dB;
targetResults=allResults(targetMask,:);
targetResults.MeetsBER=targetResults.CoopBER<=targetBER;
targetResults.MeetsOutage=targetResults.CoopOutage<=targetOutage;
targetResults.Feasible=targetResults.MeetsBER & targetResults.MeetsOutage;
writetable(targetResults,fullfile(outputFolder,'partD_target_comparison.csv'));
feasible=find(targetResults.Feasible);
hasRecommendation=~isempty(feasible);
if hasRecommendation
    % Lowest measured BER, then outage; highest effective throughput breaks
    % any remaining tie. This ranking is conditional on stated targets.
    score=[targetResults.CoopBER(feasible), ...
        targetResults.CoopOutage(feasible), ...
        -targetResults.CoopThroughput(feasible)];
    [~,order]=sortrows(score,[1 2 3]);
    winner=targetResults(feasible(order(1)),:);
else
    winner=targetResults([],:); % no design is reported as target-compliant
end

%% Part D figures
colors=lines(nPos);
figure('Color','w'); hold on;
plot(EbN0dB,C.BER_direct,'k-o','LineWidth',1.4);
for p=1:nPos
    plot(EbN0dB,reshape(BER(defaultIndex,p,:),1,nSNR),'-s', ...
        'Color',colors(p,:),'LineWidth',1.4);
end
yline(targetBER,'k:','LineWidth',1.2);
grid on; xlabel('Nominal direct-reference E_b/N_0 (dB)');
ylabel('Measured payload BER');
title('Where the 50/50 cooperative link improves reliability');
legend([{'Direct'},arrayfun(@(v)sprintf('Relay %g m',v), ...
    positions,'UniformOutput',false),{'BER target'}], ...
    'Location','northeast');
exportgraphics(gcf,fullfile(outputFolder,'partD_reliability_regions.png'), ...
    'Resolution',200);

% Separate figures keep the axis labels readable in the handwritten report.
figure('Color','w','Position',[100 100 850 500]); hold on;
for p=1:nPos
    plot(sourceFraction,reshape(BER(:,p,targetIndex),1,nSplit),'-o', ...
        'Color',colors(p,:),'LineWidth',1.5);
end
yline(targetBER,'k:'); grid on;
xticks(sourceFraction); xticklabels({'25/75','50/50','75/25'});
xlim([0.18 0.82]);
xlabel('Source/relay energy split (%)');
ylabel('Payload BER at target E_b/N_0'); title('Power split versus BER');
legend([arrayfun(@(v)sprintf('%g m',v),positions, ...
    'UniformOutput',false),{'Target'}],'Location','best');
exportgraphics(gcf,fullfile(outputFolder,'partD_power_split_BER.png'), ...
    'Resolution',200);

figure('Color','w','Position',[100 100 850 500]); hold on;
for p=1:nPos
    plot(sourceFraction,reshape(OUT(:,p,targetIndex),1,nSplit),'-s', ...
        'Color',colors(p,:),'LineWidth',1.5);
end
yline(targetOutage,'k:'); grid on;
xticks(sourceFraction); xticklabels({'25/75','50/50','75/25'});
xlim([0.18 0.82]);
xlabel('Source/relay energy split (%)');
ylabel('Outage probability at target E_b/N_0');
title('Power split versus outage');
legend([arrayfun(@(v)sprintf('%g m',v),positions, ...
    'UniformOutput',false),{'Target'}],'Location','best');
exportgraphics(gcf,fullfile(outputFolder,'partD_power_split_outage.png'), ...
    'Resolution',200);
% Remove the older two-panel figure if this script is rerun in the same folder.
oldFigure=fullfile(outputFolder,'partD_power_split.png');
if isfile(oldFigure), delete(oldFigure); end

% This small table is the key for individual points in the trade-off plot.
SourcePercent=100*targetResults.Alpha;
RelayPercent=100*(1-targetResults.Alpha);
plotKey=table(targetResults.PositionM,SourcePercent,RelayPercent, ...
    targetResults.CoopBER,targetResults.CoopOutage, ...
    targetResults.CoopThroughput, ...
    'VariableNames',{'RelayPositionM','SourcePercent','RelayPercent', ...
    'BER','Outage','UsefulBitsPerSymbolPeriod'});
writetable(plotKey,fullfile(outputFolder,'partD_tradeoff_key.csv'));

figure('Color','w','Position',[100 100 850 550]); hold on;
for p=1:nPos
    rows=targetResults.PositionM==positions(p);
    scatter(targetResults.CoopThroughput(rows), ...
        max(targetResults.CoopBER(rows),1/(nFrames*nPayload)), ...
        85,colors(p,:),'filled');
end
scatter(C.TH_direct(targetIndex),max(C.BER_direct(targetIndex), ...
    1/(nFrames*nPayload)),95,'k','filled');
yline(targetBER,'k:'); set(gca,'YScale','log'); grid on;
xlabel('Delivered payload bits per symbol period');
ylabel('Payload BER (zero shown at one-bit display floor)');
title(sprintf('Reliability versus slot efficiency at %g dB (splits in table)',targetEbN0dB));
legend([arrayfun(@(v)sprintf('Cooperative %g m',v),positions, ...
    'UniformOutput',false),{'Direct','BER target'}], ...
    'Location','best');
exportgraphics(gcf,fullfile(outputFolder,'partD_tradeoff.png'), ...
    'Resolution',200);

%% Numeric narrative for the handwritten engineering decision.
reportFile=fullfile(outputFolder,'partD_decision.txt');
fid=fopen(reportFile,'w'); assert(fid~=-1);
fprintf(fid,'PART D - EMPIRICAL ENGINEERING DECISION\n');
fprintf(fid,'Target: BER <= %.4g AND symbol-SNR outage <= %.4g at nominal %.1f dB.\n', ...
    targetBER,targetOutage,targetEbN0dB);
fprintf(fid,'Outage uses instantaneous post-combining symbol SNR < %g dB.\n', ...
    outageThresholdDB);
fprintf(fid,'Energy cap: %.0f per attempted frame; source shares %s.\n', ...
    totalEnergy,mat2str(sourceFraction));
fprintf(fid,'Cooperation reserves two 500-symbol slots; direct uses one.\n');
fprintf(fid,'Direct at target: BER %.6g, PER %.3f, outage %.3f, throughput %.3f bits/symbol-period.\n', ...
    C.BER_direct(targetIndex),C.PER_direct(targetIndex), ...
    C.OUT_direct(targetIndex),C.TH_direct(targetIndex));
fprintf(fid,'\nDefault 50/50 split: BER comparison over the %d tested Eb/N0 values.\n',nSNR);
for p=1:nPos
    b=reshape(BER(defaultIndex,p,:),1,nSNR);
    wins=EbN0dB(b<C.BER_direct);
    nonwins=EbN0dB(b>=C.BER_direct);
    fprintf(fid,'  Relay %g m: improves at %s dB; does not improve at %s dB.\n', ...
        positions(p),mat2str(wins),mat2str(nonwins));
end
fprintf(fid,'\nAt target Eb/N0 (source/relay energy, relay position, BER, outage, throughput, passes):\n');
for i=1:height(targetResults)
    row=targetResults(i,:);
    fprintf(fid,'  %g/%g, %g m: %.6g, %.3f, %.3f, %d\n', ...
        100*row.Alpha,100*(1-row.Alpha),row.PositionM, ...
        row.CoopBER,row.CoopOutage,row.CoopThroughput,row.Feasible);
end
if hasRecommendation
    fprintf(fid,'\nRECOMMENDATION UNDER THE STATED TARGET: relay %g m, ',winner.PositionM);
    fprintf(fid,'source/relay split %g/%g.\n', ...
        100*winner.Alpha,100*(1-winner.Alpha));
    fprintf(fid,'Observed BER %.6g, PER %.3f, outage %.3f, throughput %.3f bits/symbol-period.\n', ...
        winner.CoopBER,winner.CoopPER,winner.CoopOutage, ...
        winner.CoopThroughput);
else
    fprintf(fid,'\nNO CONFIGURATION MEETS BOTH TARGETS at this Eb/N0.\n');
    fprintf(fid,'Raise Eb/N0, reconsider geometry/budget or revise the target; do not claim a feasible recommendation.\n');
end
fprintf(fid,'\nLimits: %d frames per condition; PER/outage resolution %.3f.\n', ...
    nFrames,1/nFrames);
fprintf(fid,'This is a reproducible simulated comparison, not a field measurement or confidence guarantee.\n');
fclose(fid);
type(reportFile);

save(fullfile(outputFolder,'partD_data.mat'),'sourceFraction','positions', ...
    'EbN0dB','BER','PER','OUT','TH','forwardRate','allResults', ...
    'targetResults','winner','targetBER','targetOutage', ...
    'targetEbN0dB','undetectedRelay','hasRecommendation');

function bits=qpsk_demap_d(symbols)
bits=zeros(2*numel(symbols),1);
bits(1:2:end)=real(symbols(:))<0;
bits(2:2:end)=imag(symbols(:))<0;
end

function symbols=qpsk_map_d(bits)
bits=bits(:);
symbols=((1-2*bits(1:2:end))+1i*(1-2*bits(2:2:end)))/sqrt(2);
end

function valid=crc_ok_d(frame)
frame=frame(:);
valid=isequal(frame(end-7:end),crc8_bits_d(frame(1:end-8)));
end

function remainder=crc8_bits_d(bits)
% Part A CRC-8, polynomial 0x07, init 0, no reflection or final XOR.
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

%% Part B - three-link propagation model for Question 16 (base MATLAB)
% Put this script beside partA_data.mat, then run it after Part A.
% This models the channels only. Part C will decode at the relay and form
% the relay-to-destination observation using hRD and noiseRD saved here.
clear; clc; close all;
outputFolder=fileparts(mfilename('fullpath'));
partAFile=fullfile(outputFolder,'partA_data.mat');
assert(isfile(partAFile),'Run Part A first and put partA_data.mat beside this script.');
loaded=load(partAFile,'tx');
assert(isfield(loaded,'tx'),'partA_data.mat must contain tx.');
tx=loaded.tx;
assert(size(tx.sourceTx,1)==tx.symbolsPerFrame);
assert(size(tx.sourceTx,2)==tx.cfg.numberOfFrames);
assert(abs(tx.sourceEnergy+tx.relayEnergy-tx.cfg.totalEnergyPerFrame)<1e-10);

%% Geometry and propagation assumptions
cfg.seed=3128;
cfg.sourcePositionM=0;
cfg.destinationPositionM=100;
cfg.relayPositionsM=[25 50 75];
cfg.referenceDistanceM=1;
cfg.referenceLossDB=30;           % Illustrative assumed loss at 1 m
cfg.pathLossExponent=3;
cfg.EbN0dB=0:2:14;               % Nominal full-budget direct-link reference
cfg.twoBitsPerQPSKSymbol=2;
cfg.fading='Independent block Rayleigh: constant within a frame';
cfg.channelKnowledge='Part C must estimate or assume knowledge of h';
cfg.pathLossConvention='Relative to the 100 m S-D link';

nSymbols=tx.symbolsPerFrame;
nFrames=tx.cfg.numberOfFrames;
nSNR=numel(cfg.EbN0dB);
nPos=numel(cfg.relayPositionsM);
dSD=cfg.destinationPositionM-cfg.sourcePositionM;
dSR=cfg.relayPositionsM-cfg.sourcePositionM;
dRD=cfg.destinationPositionM-cfg.relayPositionsM;
assert(all(dSR>0 & dRD>0));

% Physical path-loss values for the link-budget table. The waveform model
% is normalized so that S-D power gain is 1, avoiding an unstated receiver
% noise figure or bandwidth. Apply the same 90 dB normalization to each link.
loss=@(d) cfg.referenceLossDB+10*cfg.pathLossExponent* ...
    log10(d/cfg.referenceDistanceM);
PL_SD_dB=loss(dSD);
PL_SR_dB=loss(dSR);
PL_RD_dB=loss(dRD);
gainSD=1;
gainSR=10.^((PL_SD_dB-PL_SR_dB)/10);
gainRD=10.^((PL_SD_dB-PL_RD_dB)/10);
assert(max(abs(gainSR-(dSD./dSR).^cfg.pathLossExponent))<1e-11);
assert(max(abs(gainRD-(dSD./dRD).^cfg.pathLossExponent))<1e-11);

linkBudget=table(cfg.relayPositionsM(:),dSR(:),dRD(:), ...
    repmat(dSD,nPos,1),repmat(PL_SD_dB,nPos,1), ...
    PL_SR_dB(:),PL_RD_dB(:),gainSR(:),gainRD(:), ...
    'VariableNames',{'RelayPositionM','DistanceSR_M','DistanceRD_M', ...
    'DistanceSD_M','PathLossSD_dB','PathLossSR_dB', ...
    'PathLossRD_dB','RelativeGainSR','RelativeGainRD'});
writetable(linkBudget,fullfile(outputFolder,'partB_link_budget.csv'));
disp(linkBudget);

%% Noise and independent fading for every frame and Eb/N0 setting
% A full-budget one-slot direct transmission spends Etotal/nSymbols per
% unit-energy QPSK symbol. Thus Eb_ref=Es_ref/2 and N0=Eb_ref/(Eb/N0).
% Complex noise has E{|n|^2}=N0, with N0/2 in each real dimension.
EsReference=tx.cfg.totalEnergyPerFrame/nSymbols;
N0=EsReference./(cfg.twoBitsPerQPSKSymbol*10.^(cfg.EbN0dB/10));
directPower=EsReference; % One-slot direct comparison, same frame energy
previousRng=rng;
restoreRng=onCleanup(@()rng(previousRng)); %#ok<NASGU>
rng(cfg.seed,'twister');

% For fair direct/cooperative comparisons, S-D h and noise are reused.
% S-R, R-D and each receiver's AWGN are otherwise independent.
hSD=zeros(nFrames,nSNR);
hSR=zeros(nFrames,nSNR,nPos);
hRD=zeros(nFrames,nSNR,nPos);
rxSD=zeros(nSymbols,nFrames,nSNR);
rxDirect=zeros(nSymbols,nFrames,nSNR);
rxSR=zeros(nSymbols,nFrames,nSNR,nPos);
noiseRD=zeros(nSymbols,nFrames,nSNR,nPos);
snrSD=zeros(nFrames,nSNR);
snrDirect=zeros(nFrames,nSNR);
snrSR=zeros(nFrames,nSNR,nPos);
snrRD=zeros(nFrames,nSNR,nPos);

for k=1:nSNR
    noiseScale=sqrt(N0(k)/2);
    for f=1:nFrames
        hSD(f,k)=(randn+1i*randn)/sqrt(2);
        nSD=noiseScale*(randn(nSymbols,1)+1i*randn(nSymbols,1));
        rxSD(:,f,k)=sqrt(gainSD)*hSD(f,k)*tx.sourceTx(:,f)+nSD;
        rxDirect(:,f,k)=sqrt(gainSD*directPower)* ...
            hSD(f,k)*tx.sourceSymbols(:,f)+nSD;
        snrSD(f,k)=tx.sourcePower*gainSD*abs(hSD(f,k))^2/N0(k);
        snrDirect(f,k)=directPower*gainSD*abs(hSD(f,k))^2/N0(k);
        for p=1:nPos
            hSR(f,k,p)=(randn+1i*randn)/sqrt(2);
            hRD(f,k,p)=(randn+1i*randn)/sqrt(2);
            nSR=noiseScale*(randn(nSymbols,1)+1i*randn(nSymbols,1));
            noiseRD(:,f,k,p)=noiseScale* ...
                (randn(nSymbols,1)+1i*randn(nSymbols,1));
            rxSR(:,f,k,p)=sqrt(gainSR(p))* ...
                hSR(f,k,p)*tx.sourceTx(:,f)+nSR;
            snrSR(f,k,p)=tx.sourcePower*gainSR(p)* ...
                abs(hSR(f,k,p))^2/N0(k);
            snrRD(f,k,p)=tx.relayPower*gainRD(p)* ...
                abs(hRD(f,k,p))^2/N0(k);
        end
    end
end

% Part C, after CRC-controlled decode-and-forward, should compute:
% yRD(:,f,k,p)=sqrt(gainRD(p))*hRD(f,k,p)*relayTx(:,f) ...
%                 +noiseRD(:,f,k,p);
% relayTx already contains sqrt(tx.relayPower) and is zero on CRC failure.
% On relay failure, omit the R-D branch rather than combining pure noise.

assert(all(isfinite(rxSR(:))) && all(isfinite(rxSD(:))));
assert(all(abs(N0-EsReference./(2*10.^(cfg.EbN0dB/10)))<1e-12));

%% Part B evidence plots
figure('Color','w');
plot(cfg.relayPositionsM,PL_SR_dB,'o-','LineWidth',1.6); hold on;
plot(cfg.relayPositionsM,PL_RD_dB,'s-','LineWidth',1.6);
yline(PL_SD_dB,'k:','LineWidth',1.4);
grid on; xlabel('Relay position from source (m)'); ylabel('Path loss (dB)');
title('Three-link distance-based path loss');
legend('Source to relay','Relay to destination','Source to destination', ...
    'Location','best');
exportgraphics(gcf,fullfile(outputFolder,'partB_path_loss.png'),'Resolution',200);

kExample=find(cfg.EbN0dB==8,1);
figure('Color','w');
plot(1:nFrames,abs(hSD(:,kExample)).^2,'k-','LineWidth',1.1); hold on;
plot(1:nFrames,abs(hSR(:,kExample,2)).^2,'b-','LineWidth',1.1);
plot(1:nFrames,abs(hRD(:,kExample,2)).^2,'r-','LineWidth',1.1);
grid on; xlabel('Frame number'); ylabel('Rayleigh power gain |h|^2');
title('Independent block fading, relay at 50 m');
legend('S-D','S-R','R-D','Location','best');
exportgraphics(gcf,fullfile(outputFolder,'partB_fading.png'),'Resolution',200);

figure('Color','w');
plot(cfg.EbN0dB,10*log10(mean(snrDirect,1)),'k-o','LineWidth',1.5);
hold on;
for p=1:nPos
    plot(cfg.EbN0dB,10*log10(mean(snrSR(:,:,p),1)),'-o', ...
        'LineWidth',1.2);
    plot(cfg.EbN0dB,10*log10(mean(snrRD(:,:,p),1)),'--s', ...
        'LineWidth',1.2);
end
grid on; xlabel('Nominal direct-link E_b/N_0 (dB)');
ylabel('Mean instantaneous link symbol SNR (dB)');
title('Channel conditions before relay/destination decoding');
legend('Direct S-D','S-R 25 m','R-D 25 m','S-R 50 m', ...
    'R-D 50 m','S-R 75 m','R-D 75 m','Location','best');
exportgraphics(gcf,fullfile(outputFolder,'partB_link_snr.png'),'Resolution',200);

%% Save sufficient observations for CRC forwarding and MRC in Part C
save(fullfile(outputFolder,'partB_data.mat'),'cfg','N0','EsReference', ...
    'linkBudget','gainSD','gainSR','gainRD','hSD','hSR','hRD', ...
    'rxSD','rxDirect','rxSR','noiseRD','snrSD','snrDirect', ...
    'snrSR','snrRD','-v7');
fid=fopen(fullfile(outputFolder,'partB_validation.txt'),'w');
assert(fid~=-1);
fprintf(fid,'PART B CHANNEL VALIDATION\n');
fprintf(fid,'Source 0 m, destination 100 m; relays 25, 50, 75 m.\n');
fprintf(fid,'Log-distance exponent 3; assumed reference loss 30 dB at 1 m.\n');
fprintf(fid,'Direct path loss %.2f dB, normalized direct power gain 1.\n',PL_SD_dB);
fprintf(fid,'Independent block Rayleigh links and AWGN; seed %d.\n',cfg.seed);
fprintf(fid,'Nominal direct-reference Eb/N0: %s dB.\n',mat2str(cfg.EbN0dB));
fprintf(fid,'Noise N0 = (full-budget symbol energy / 2) / nominal Eb/N0.\n');
fprintf(fid,'Direct baseline spends full frame energy in one time slot.\n');
fprintf(fid,'Relay-phase received symbols await Part C CRC forwarding.\n');
fprintf(fid,'No BER, PER, outage, throughput or optimum placement is claimed here.\n');
fclose(fid);
type(fullfile(outputFolder,'partB_validation.txt'));

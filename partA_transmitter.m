%% Run Part A demonstration and export evidence (base MATLAB only).
% Keep all companion .m files in this folder.
clear; clc; close all;
outputFolder=fileparts(mfilename('fullpath'));
addpath(outputFolder);
tx=transmitter();

%% Test 1: exact noiseless recovery of every payload and complete frame
idealRelay=relay_decode_forward(tx.sourceTx,tx);
payloadErrors=nnz(idealRelay.payloadBits~=tx.payloadBits);
assert(payloadErrors==0);
assert(isequal(idealRelay.detectedFrames,tx.sourceFrames));
assert(all(idealRelay.forward));
assert(isequal(qpsk_codec(idealRelay.relayTx,'demap'),tx.sourceFrames));
assert(abs(mean(abs(tx.sourceSymbols(:)).^2)-1)<1e-12);

%% Test 2: separate controlled CRC rejection, not channel performance
corruptedInput=tx.sourceTx;
% Flip only the first payload bit (symbol 9, real component).
corruptedInput(9,1)=-real(corruptedInput(9,1))+1j*imag(corruptedInput(9,1));
crcTestRelay=relay_decode_forward(corruptedInput,tx);
assert(~crcTestRelay.forward(1) && all(crcTestRelay.forward(2:end)));
assert(all(crcTestRelay.relayTx(:,1)==0));
assert(nnz(crcTestRelay.detectedFrames~=tx.sourceFrames)==1);

%% Per-frame budget and nonoverlapping time slots
energyUsed=sum(abs(tx.sourceTx).^2,1)+sum(abs(idealRelay.relayTx).^2,1);
assert(all(abs(energyUsed-tx.cfg.totalEnergyPerFrame)<1e-9));
assert(all(tx.sourceTwoPhase.*idealRelay.relayTwoPhase==0,'all'));
crcTestEnergy=sum(abs(tx.sourceTx).^2,1)+sum(abs(crcTestRelay.relayTx).^2,1);
assert(all(crcTestEnergy<=tx.cfg.totalEnergyPerFrame+1e-9));

%% Evidence: time, energy, constellation and payload
figure('Color','w');
slotTime=[0 tx.symbolsPerFrame 2*tx.symbolsPerFrame];
stairs(slotTime,[tx.sourcePower 0 0],'b-','LineWidth',2); hold on;
stairs(slotTime,[0 tx.relayPower 0],'r--','LineWidth',2);
xticks(slotTime); xlim([0 slotTime(end)]);
ylim([0 1.2*max(tx.sourcePower,tx.relayPower)]);
xlabel('Time in symbol periods (phase 1 then phase 2)');
ylabel('Normalized transmit power'); legend('Source','Relay');
title('Two-phase transmission: a CRC-valid frame'); grid on;
exportgraphics(gcf,fullfile(outputFolder,'partA_timing.png'),'Resolution',200);

figure('Color','w');
bar([tx.sourceEnergy tx.relayEnergy;tx.sourceEnergy 0],'stacked');
xticks([1 2]); xticklabels({'CRC passes','CRC fails'});
ylabel('Energy per attempted frame (normalized units)');
legend('Source','Relay'); title('Fixed energy budget; unused reserve on failure'); grid on;
exportgraphics(gcf,fullfile(outputFolder,'partA_energy.png'),'Resolution',200);

figure('Color','w');
plot(real(tx.sourceSymbols(:)),imag(tx.sourceSymbols(:)),'o');
xlabel('In-phase'); ylabel('Quadrature'); axis equal;
xlim([-1 1]); ylim([-1 1]); grid on; title('Unit-energy QPSK mapping');
exportgraphics(gcf,fullfile(outputFolder,'partA_constellation.png'),'Resolution',200);

figure('Color','w');
stairs(1:32,tx.payloadBits(1:32,1),'LineWidth',1.5);
ylim([-0.2 1.2]); yticks([0 1]); xlim([1 32]); grid on;
xlabel('Payload bit index'); ylabel('Bit value');
title('First 32 payload bits of frame ID 1 (seed 16)');
exportgraphics(gcf,fullfile(outputFolder,'partA_payload.png'),'Resolution',200);

%% Configurable energy splits, not performance optimization
alpha=[0.25;0.50;0.75]; budget=tx.cfg.totalEnergyPerFrame;
SourceEnergy=alpha*budget; RelayEnergy=(1-alpha)*budget;
SourcePower=SourceEnergy/tx.symbolsPerFrame;
RelayPower=RelayEnergy/tx.symbolsPerFrame;
powerTable=table(alpha,SourceEnergy,RelayEnergy,SourcePower,RelayPower);
writetable(powerTable,fullfile(outputFolder,'partA_power_allocation.csv'));
for k=1:numel(alpha)
    variant=transmitter(struct('sourceEnergyFraction',alpha(k)));
    forwarded=relay_decode_forward(variant.sourceTx,variant);
    actual=sum(abs(variant.sourceTx).^2,1)+sum(abs(forwarded.relayTx).^2,1);
    assert(all(abs(actual-budget)<1e-9));
end

testDescription='Ideal baseline and separate one-bit CRC test; no channel model or measured link performance.';
save(fullfile(outputFolder,'partA_data.mat'),'tx','idealRelay','crcTestRelay', ...
    'testDescription','powerTable','payloadErrors');
fid=fopen(fullfile(outputFolder,'partA_validation.txt'),'w'); assert(fid~=-1);
fprintf(fid,'PART A FUNCTIONAL VALIDATION\n');
fprintf(fid,'Frame: 16 ID + 976 payload + 8 CRC = 1000 bits; 500 QPSK symbols.\n');
fprintf(fid,'Ideal test: %d payload errors across %d bits; 100/100 CRC pass.\n', ...
    payloadErrors,numel(tx.payloadBits));
fprintf(fid,'Entire frames and forwarded frames recovered exactly.\n');
fprintf(fid,'Separate CRC test: one injected bit error; frame 1 rejected; 99 other frames forwarded.\n');
fprintf(fid,'Default energy: 250 source + 250 relay = 500 per attempted frame.\n');
fprintf(fid,'CRC rejection uses 250 units; 250 units remain unused.\n');
fprintf(fid,'25/75, 50/50 and 75/25 energy splits passed per-frame budget checks.\n');
fprintf(fid,'No channel, destination combining or link-performance metrics are modeled.\n');
fclose(fid);
type(fullfile(outputFolder,'partA_validation.txt'));
disp(powerTable);

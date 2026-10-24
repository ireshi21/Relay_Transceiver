%% Run Part A demonstration and export evidence (base MATLAB only).
% Self-contained script: all helper functions are included below.
clear; clc; close all;
outputFolder=fileparts(mfilename('fullpath'));
% No companion MATLAB files are required.
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


function tx = transmitter(cfg)
% Part A source module. Columns are frames of complex baseband symbols.
% Format: 16-bit ID (MSB first) | payload | CRC-8 of ID plus payload.
if nargin==0, cfg=struct(); end
defaults=struct('seed',16,'numberOfFrames',100,'payloadBitsPerFrame',976, ...
    'totalEnergyPerFrame',500,'sourceEnergyFraction',0.5);
names=fieldnames(defaults);
for k=1:numel(names)
    if ~isfield(cfg,names{k}), cfg.(names{k})=defaults.(names{k}); end
end
validateattributes(cfg.numberOfFrames,{'numeric'},{'scalar','integer','>=',1,'<=',65535});
validateattributes(cfg.payloadBitsPerFrame,{'numeric'},{'scalar','integer','positive'});
validateattributes(cfg.sourceEnergyFraction,{'numeric'},{'scalar','>',0,'<',1});
validateattributes(cfg.totalEnergyPerFrame,{'numeric'},{'scalar','real','finite','positive'});
assert(mod(cfg.payloadBitsPerFrame,2)==0,'Payload length must be even.');
previousRng=rng; restoreRng=onCleanup(@()rng(previousRng)); %#ok<NASGU>
rng(cfg.seed,'twister');
tx.cfg=cfg; tx.frameIdBits=16; tx.crcBits=8;
tx.frameBits=16+cfg.payloadBitsPerFrame+8;
tx.symbolsPerFrame=tx.frameBits/2;
tx.payloadRows=17:(16+cfg.payloadBitsPerFrame);
tx.payloadBits=randi([0 1],cfg.payloadBitsPerFrame,cfg.numberOfFrames);
tx.frameIds=1:cfg.numberOfFrames;
tx.sourceFrames=zeros(tx.frameBits,cfg.numberOfFrames);
for k=1:cfg.numberOfFrames
    id=double(bitget(uint16(k),16:-1:1)).';
    body=[id;tx.payloadBits(:,k)];
    tx.sourceFrames(:,k)=[body;crc8_bits(body)];
end
tx.sourceSymbols=qpsk_codec(tx.sourceFrames,'map');
tx.sourceEnergy=cfg.sourceEnergyFraction*cfg.totalEnergyPerFrame;
tx.relayEnergy=(1-cfg.sourceEnergyFraction)*cfg.totalEnergyPerFrame;
% Normalized symbol duration Ts=1; energy=sum(abs(symbols).^2).
tx.sourcePower=tx.sourceEnergy/tx.symbolsPerFrame;
tx.relayPower=tx.relayEnergy/tx.symbolsPerFrame;
tx.sourceTx=sqrt(tx.sourcePower)*tx.sourceSymbols;
tx.sourceTwoPhase=[tx.sourceTx;zeros(size(tx.sourceTx))];
end


function relay=relay_decode_forward(equalizedSymbols,tx)
% Part A relay strategy. Input: coherent equalized symbol columns.
% Source amplitude scaling may remain; phase rotation must be corrected.
% Output is already scaled by relay power. No propagation model is included.
assert(isequal(size(equalizedSymbols),size(tx.sourceTx)));
assert(all(isfinite(equalizedSymbols),'all'));
relay.detectedFrames=qpsk_codec(equalizedSymbols,'demap');
relay.forward=false(1,size(equalizedSymbols,2));
relay.relayTx=zeros(size(equalizedSymbols));
for k=1:size(equalizedSymbols,2)
    frame=relay.detectedFrames(:,k);
    relay.forward(k)=isequal(crc8_bits(frame(1:end-8)),frame(end-7:end));
    if relay.forward(k)
        relay.relayTx(:,k)=sqrt(tx.relayPower)*qpsk_codec(frame,'map');
    end
end
relay.payloadBits=relay.detectedFrames(tx.payloadRows,:);
relay.relayTwoPhase=[zeros(size(relay.relayTx));relay.relayTx];
end


function out=qpsk_codec(in,mode)
% Shared mapping: 00 -> (1+j)/sqrt(2); 01 -> (1-j)/sqrt(2);
% 10 -> (-1+j)/sqrt(2); 11 -> (-1-j)/sqrt(2).
switch mode
    case 'map'
        assert(mod(size(in,1),2)==0 && all(in==0 | in==1,'all'));
        out=((1-2*in(1:2:end,:))+1j*(1-2*in(2:2:end,:)))/sqrt(2);
    case 'demap'
        out=zeros(2*size(in,1),size(in,2));
        out(1:2:end,:)=real(in)<0;
        out(2:2:end,:)=imag(in)<0;
    otherwise
        error('Mode must be map or demap.');
end
end


function remainder=crc8_bits(bits)
% MSB-first division, g(x)=x^8+x^2+x+1 (0x07).
% Initial state zero, no reflection, no final XOR; append eight zeros.
bits=double(bits(:));
assert(all(bits==0 | bits==1),'CRC input must be binary.');
work=[bits;zeros(8,1)]; generator=[1 0 0 0 0 0 1 1 1]';
for i=1:numel(bits)
    if work(i)==1
        work(i:i+8)=mod(work(i:i+8)+generator,2);
    end
end
remainder=work(end-7:end);
end


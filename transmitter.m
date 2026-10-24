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

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

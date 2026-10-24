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

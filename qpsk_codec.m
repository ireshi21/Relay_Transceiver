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

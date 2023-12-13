%convert bits to QPSK symbols
function info_symbols=bits2QPSK(info_bits)
if rem(length(info_bits),2)~=0
    info_bits=[info_bits 0];
end
L=length(info_bits)/2;
info_symbols=zeros(L,1);
for k=1:L;
    info_symbols(k)=2*(info_bits(2*k-1)-0.5)+j*2*(info_bits(2*k)-0.5);
end
    
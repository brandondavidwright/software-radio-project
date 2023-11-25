%
% Analysis filter bank for CMT (AFT-CMT) 
%
% Inputs:
%   y: the signal to be analyzed
%   L: Oversampling factor (size of ifft = 2L)
%   h: Prototype filter
%
% Output:
%   Y: Analyzed signal 
%
function Y=AFB_CMT(y,L,h);
w2L=exp(-j*pi/L);
outputE=[];
for k=1:2*L
    temp=h(k:2*L:end);
    for m=1:size(temp)
        temp(m)=(-1)^(m+1)*temp(m);
    end
    Ek=expander(temp,2);
    yk=[zeros(1,k-1) (w2L^(-1/2*(k-1))*y).' zeros(1,2*L-k+1)];
    outputE=[outputE;conv(Ek,yk(1:L:end))]; 
end

Y=2*L*ifft(outputE,2*L);

%
% Analysis filter bank (AFT) 
%
% Inputs:
%   x: the signal to be analyzed
%   L: Oversampling factor (size of ifft)
%   h: Prototype filter
%
% Output:
%   S: Analyzed/extracted subcarrier signals/symbols 
%
function S=AFB(x,L,h);
X=[];
x=[x;zeros(L-mod(length(x),L),1)];
for k=1:L
    temp=[zeros(k-1,1);x];
    temp=temp(1:L:end-k);
    X=[X;(conv(h(k:L:end),temp)).'];
end
S=L*ifft(X,L);


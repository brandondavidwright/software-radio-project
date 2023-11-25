%
% Synthesis filter bank (SFB) 
%
% Inputs:
%   S: Matrix of data symbols to be synthesized to a multicarrier signal
%   L: Oversampling factor (size of ifft)
%   h: Prototype filter
%
% Output:
%   x: Synthesized output 
%
function x=AFB(S,L,h);
x=[];
A=L*ifft(S,L);
for k=1:L
    x=[x;conv(h(k:L:end),A(k,:))];
end
x=reshape(x,[],1);

%
% Synthesis filter bank for CMT (SFB-CMT) 
%
% Inputs:
%   A: Matrix of data symbols scaled with j^k
%   L: Oversampling factor (size of ifft = 2L)
%   h: Prototype filter
%
% Output:
%   x: Synthesized output 
%
function x=SFB_CMT(A,L,h);

IDFToutput=2*L*ifft(A,2*L);
outputE=[];
w2L=exp(-j*pi/L);
for k=1:2*L
    temp=h(k:2*L:end);
    for m=1:size(temp)
        temp(m)=(-1)^(m+1)*temp(m);
    end
    Ek=expander(temp,2);
    outputE=[outputE;conv(Ek,IDFToutput(k,:))]; 
end

x=zeros(1,size(outputE,2)*L+2*L-1);
for k=1:2*L
    temp2=[zeros(1,k-1) w2L^(-1/2*(k-1))*expander(outputE(k,:),L) zeros(1,2*L-k)];
    x=x+temp2;
end

x=x.';

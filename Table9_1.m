phic%
% This program can be used to generate the results of Table 9.1 of the
% text.
%
%%%%%%%%%%%%%%%%%%%%%%
%     PARAMETERS     %
%%%%%%%%%%%%%%%%%%%%%%
clear all
Tb=0.0001;      % Symbol/baud period 
L=5;          % Number of samples per symbol period
Ts=Tb/L;        % Sampling period
fs=1/Ts;
fc=100000;      % Carrier frequency at the transmitter
Dfc=10;      % Carrier frequency offset
phic=0;        % Carrier phase offset
alpha=1;      % Roll-off factor for the (square-root) raised cosine filters 
K=8;            % Half of the length of the square-root raised cosine filters in Tb. 
sigma_v=0;      % Standard deviation of channel noise

%%%%%%%%%%%%%%%%%%%%%%
%       SOURCE       %
%%%%%%%%%%%%%%%%%%%%%%
N=600000;
b=sign(randn(N,1)); 
%%%%%%%%%%%%%%%%%%%%%%
%   Transmit filter  %
%%%%%%%%%%%%%%%%%%%%%%
pT=sr_cos_p(K*L,L,alpha);


%%%%%%%%%%%%%%%%%%%%%%
%      4-QAM/QPSK    %
%%%%%%%%%%%%%%%%%%%%%%
s=b(1:2:end)+i*b(2:2:end);
x=conv(expander(s,L),pT);
ym4(1)=mean(real(x.^4))/mean(abs(x).^4);




%%%%%%%%%%%%%%%%%%%%%%
%      16-QAM        %
%%%%%%%%%%%%%%%%%%%%%%
s=zeros(length(b)/4,1);
for k=1:length(b)/4;                        
    kk=(k-1)*4;
    s(k)=2*b(kk+1)+b(kk+2)+j*(2*b(kk+3)+b(kk+4));
end
x=conv(expander(s,L),pT);
ym4(2)=mean(real(x.^4))/mean(abs(x).^4);


%%%%%%%%%%%%%%%%%%%%%%
%      64-QAM        %
%%%%%%%%%%%%%%%%%%%%%%
s=zeros(length(b)/6,1);
for k=1:length(b)/6;                        
    kk=(k-1)*6;
    s(k)=4*b(kk+1)+2*b(kk+2)+b(kk+3)+j*(4*b(kk+4)+2*b(kk+5)+b(kk+6));
end
x=conv(expander(s,L),pT);
ym4(3)=mean(real(x.^4))/mean(abs(x).^4);




%%%%%%%%%%%%%%%%%%%%%%
%      256-QAM       %
%%%%%%%%%%%%%%%%%%%%%%
s=zeros(length(b)/8,1);
for k=1:length(b)/8;                        
    kk=(k-1)*8;
    s(k)=8*b(kk+1)+4*b(kk+2)+2*b(kk+3)+b(kk+4)+j*(8*b(kk+5)+4*b(kk+6)+2*b(kk+7)+b(kk+8));
end
x=conv(expander(s,L),pT);
ym4(4)=mean(real(x.^4))/mean(abs(x).^4);

alpha
ym4

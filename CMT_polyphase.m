%
%   CMT Transceiver 
%
clear all; close all;

Nactive=96; % Number of active subcarriers
Nsymbols=100; Tb=0.0001; L=128; Ts=Tb/L;
fc=50000; Dfc=0; phic=0; K=8; sigmav=0; c=1;

%%% SOURCE %%%
b=sign(randn(2*Nactive*Nsymbols,1));

%%% CODER %%%
s=b; % Convert bits to PAM symbols

%%% Multiplexing
A=reshape(s,Nactive,length(s)/Nactive);
A=[A;zeros(2*L-Nactive,size(A,2))];


%%% Prototype filter design %%%
K=8;M=2*L;N=K*M;alpha=1;gamma=1;
h=sr_Nyquist_p(N,M,alpha,gamma);
h=[h;zeros(2*L-1,1)]; % Make the filter length a multiple of 2L.

%%% CMT polyphase %%%
A1=[];
for k=1:size(A,1)
    A1=[A1;j^(k-1)*(A(k,:))];
end
x=SFB_CMT(A1,L,h);

%%% Channel %%%
c=[1 0 1-j 0 0.5];
y=conv(c,x);

%%% CMT receiver %%%
S=[];
Y = AFB_CMT(y,L,h);
for k=1:2*L
    S = [S;real((-j)^(k-1)*Y(k,:))];
end


sout=reshape(S(1:Nactive,:),[],1);
sout=sout(2*K*Nactive+1:end);
sout=sout(1:length(s));

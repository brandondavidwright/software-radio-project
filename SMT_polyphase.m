%
%   SMT Transceiver 
%
clear all; close all;

Nactive=48; % Number of active subcarriers
Nsymbols=100; Tb=0.0001; L=64; Ts=Tb/L;
fc=50000; Dfc=0; phic=0; K=8; sigmav=0; c=1;

%%% SOURCE %%%
b=sign(randn(2*Nactive*Nsymbols,1));

%%% CODER %%%
s=b(1:2:end)+j*b(2:2:end); % Convert bits to QPSK symbols

%%% Multiplexing
A=reshape(s,Nactive,length(s)/Nactive);
A=[A;zeros(L-Nactive,size(A,2))];

%%% Prototype filter design %%% 
K=8;M=L;N=K*M;alpha=1;gamma=1;
h=sr_Nyquist_p(N,M,alpha,gamma);
h=[h;zeros(L-1,1)]; % Make the filter length a multiple of L.  

%%% SMT transmitter %%%
A1=[];A2=[];
for k=1:size(A,1)
    A1=[A1;j^(k-1)*real(A(k,:))];
    A2=[A2;j^k*imag(A(k,:))];
end
x=[SFB(A1,L,h);zeros(L/2,1)]+[zeros(L/2,1);SFB(A2,L,h)];

x=x.*exp(j*2*pi*Dfc/L*[0:length(x)-1]');

%%% Channel %%%
c=[1];
y=conv(c,x);

%%% SMT receiver %%%
y_I = AFB(y,L,h);
y_Q = AFB([y(L/2+1:end);zeros(L/2,1)],L,h);

s_I=[];
s_Q=[];
for k=1:L
    s_I = [s_I;real((-j)^(k-1)*y_I(k,:))];
    s_Q = [s_Q;imag((-j)^(k-1)*y_Q(k,:))];
end

sout=reshape((s_I(1:Nactive,:)+j*s_Q(1:Nactive,:)),[],1);
sout=sout(K*Nactive+1:end); % removing the zeros added to the front
sout=sout(1:length(s));
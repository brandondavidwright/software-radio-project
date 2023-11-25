%
% OFDM transceiver
%
clear all, close all
%%%%%%%%%%%%%%%%
%  Parameters  %
%%%%%%%%%%%%%%%%
N=64;           % length of FFT/IFFT
Ncp=16;           % length of CP
Nactive=48;     % Number of active subcarriers
Nsymbols=50;  % Number of OFDM symbols, to transmit
Tb=0.0001;      % Symbol/baud period 
L=40;         % Number of samples per symbol period/oversampling factor
Ts=Tb/L;      % Sampling period
fc=50000;       % Carrier frequency at the transmitter
Dfc=0;          % Carrier frequency offset
phic=0;         % Carrier phase offset 
K=6;            % Length of transmit and receive filters in Tb. 
sigmav=0.001;       % Standard deviation of channel noise
c=1;   % Channel impulse response
c=[0.5; zeros(180,1); 1; zeros(88,1); 0.3];
%%%%%%%%%%
% Source %
%%%%%%%%%%
b=sign(randn(2*Nactive*Nsymbols,1));
%%%%%%%%%
% coder %
%%%%%%%%%
sfreq=b(1:2:end)+i*b(2:2:end);   % convert bits to QPSK symbols
                                 % These are in the frequency domain.
stime=OFDMmod(sfreq,N,Nactive,Ncp);   % This functions forms the OFDM symbols,
                                    % performs subcarriers modulation (IFFT),
                                    % and add cyclic prefix.
%%%%%%%%%%%%%%%%%%%%%%
% TRANSMIT FILTERING %
%%%%%%%%%%%%%%%%%%%%%%
alpha=0.25; gamma=1;
pT=sr_Nyquist_p(K*L,L,alpha,gamma);
xbbT=conv(expander(stime,L),pT);

%%%%%%%%%%%%%%%%%%%%%%
%     MODULATION     %
%%%%%%%%%%%%%%%%%%%%%%
t=[0:length(xbbT)-1]'*Ts;      % Set the time indices
xT=real(exp(i*2*pi*fc*t).*xbbT);

%%%%%%%%%%%%%%%%%%%%%%
%      CHANNEL       %
%%%%%%%%%%%%%%%%%%%%%%
xR=conv(c,xT);
xR=xR+sigmav*randn(size(xR)); % Received signal
xR=[zeros(40*80*3,1); xR];
%%%%%%%%%%%%%%%%%%%%%%
%    DEMODULATION    %
%%%%%%%%%%%%%%%%%%%%%%
t=[0:length(xR)-1]'*Ts;       % Set the time indices
xbbR=2*exp(-i*(2*pi*(fc+Dfc)*t-phic)).*xR;

%%%%%%%%%%%%%%%%%%%%%%
% RECEIVE FILTERING  %
%%%%%%%%%%%%%%%%%%%%%%
pR=pT;    
y=conv(xbbR,pR);
y=y(1:L:end);                   % Decimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Search for the best timing phase %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eta=0;
for n=1:length(y)-N-Ncp
    y1=y(n:n+Ncp-1); y2=y(n+N:n+N+Ncp-1);
    eta(n)=abs(y1'*y2)-0.5*abs(y1'*y1+y2'*y2);
end
eta=eta(1:floor(end/(N+Ncp))*(N+Ncp));
eta=reshape(eta,N+Ncp,length(eta)/(N+Ncp));
etaM=eta;
eta=mean(eta');
[etamax,k]=max(eta);k=k-K;
if k>64 k=k-64; end
if (k-1)<0 k=1; end 
figure(2),axes('position',[0.25 0.25 0.5 0.5])
plot(eta), xlabel('n'), ylabel('avg\{\eta[n]\}')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove CP, apply DFT and FEQ %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%k=1;
n=Ncp+k;
w=ones(Nactive,1); mu=0.1;
m=1; shat=[];
p=conv(pT,pR); 
c=c.*exp(-j*2*pi*[0:length(c)-1]'*Ts*fc); %c=[zeros((k-1)*L,1); c];
cBB=conv(c,p); cBB=cBB(1:L:end); cBB=[cBB; zeros(N-length(cBB),1)];
cBB=[cBB(k+1:end); zeros(N-length(cBB),1); cBB(1:k)];
CBB=fft(cBB);
while n<length(y)-N-k
    A1=fft(y(n+1:n+N));
    A1=A1./CBB;         % Equalization
    A1=[A1(2:Nactive/2+1); A1(end-Nactive/2+1:end)];
    shat=[shat; A1];
    n=n+Ncp+N;
    m=m+Nactive;
    figure(1),plot(conj(w).*A1,'.'),axis([-1.5 1.5 -1.5 1.5]),pause(0.1)
end
figure(3),plot(etaM,'k')

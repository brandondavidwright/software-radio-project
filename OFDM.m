%
% OFDM transceiver
%
clear all, close all
%%%%%%%%%%%%%%%%
%  Parameters  %
%%%%%%%%%%%%%%%%
N=64;           % length of FFT/IFFT
Ncp=16;           % length of CP
Nactive=52;     % Number of active subcarriers
Nsymbols=50;  % Number of OFDM symbols, to transmit
Tb=0.0001;      % Symbol/baud period 
L=40;         % Number of samples per symbol period/oversampling factor
Ts=Tb/L;      % Sampling period
fc=50000;       % Carrier frequency at the transmitter
Dfc=0;          % Carrier frequency offset
phic=0;         % Carrier phase offset
K=8;            % Length of transmit and receive filters in Tb. 
sigmav=0.0;       % Standard deviation of channel noise
c=1;   % Channel impulse response
%c=[1; zeros(137,1); 0.5; zeros(88,1); 0.8];
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
stime=stime+sigmav*randn(size(stime));
%%%%%%%%%%%%%%%%%%%%%%
% TRANSMIT FILTERING %
%%%%%%%%%%%%%%%%%%%%%%
pT=firpm(K*L,[0 26 38 32*L]/(32*L),[1 1 0 0],[1 10])'; 
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove CP, apply DFT and FEQ %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=Ncp;
w=ones(Nactive,1); % mu=0.1;
m=1; shat=[];
p=conv(pT,pR); c=c.*exp(-j*2*pi*[0:length(c)-1]'*Ts*fc); cBB=conv(c,p);
cBB=cBB(1:L:end); 
CBB=fft(cBB,N); 
while n<length(y)-N
    A1=fft(y(n+1:n+N));
    A1=A1./CBB;         % Equalization
    A1=[A1(2:Nactive/2+1); A1(end-Nactive/2+1:end)];

    shat=[shat; A1];
    n=n+Ncp+N;
    m=m+Nactive;
    figure(1),plot(conj(w).*A1,'.'),axis([-1.5 1.5 -1.5 1.5]),pause(0.1)
end


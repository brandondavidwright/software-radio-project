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
Nsymbols=5000;  % Number of OFDM symbols, to transmit
Tb=0.0001;      % Symbol/baud period 
L=40;         % Number of samples per symbol period/oversampling factor
Ts=Tb/L;      % Sampling period
fc=50000;       % Carrier frequency at the transmitter
Dfc=0;          % Carrier frequency offset
phic=0;         % Carrier phase offset
K=6;            % Length of transmit and receive filters in Tb. 
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
stime=OFDMmod(sfreq,N,Nactive,Ncp,1);   % This functions forms the OFDM symbols,
                                    % performs subcarriers modulation (IFFT),
                                    % and add cyclic prefix.
stime=stime+0.01*randn(size(stime));
%%%%%%%%%%%%%%%%%%%%%%
% TRANSMIT FILTERING %
%%%%%%%%%%%%%%%%%%%%%%
% pT=firpm(K*L,[0 24/32/L 40/32/L 1],[1 1 0 0],[1 10])'; 
alpha=0.22; gamma=1;
%pT=sr_Nyquist_p(K*L,L,alpha,gamma);
%pT=nqst(K*L,L,(1+alpha)/2/L);
pT=firpm(K*L-1,2*[0 26 38 32*L]/64/L,[1 1 0 0],[1 10])'; 
xbbT=conv(expander(stime,L),pT);

XX=0;
K=0;
N=2^12;
for k=1:N:length(xbbT)-N
    [X,F]=spec_analysis(xbbT(k:k+N-1).*hamming(N),1/Ts);
    XX=XX+X.^2;
    K=K+1;
end
XX=XX/max(XX);
figure,axes('position',[0.25 0.25 0.5 0.5])
plot(F,10*log10(XX))
xlabel('Frequency')
ylabel('Amplitude, dB')
axis([-2e5 2e5 -60 0])

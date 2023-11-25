close all; clear;
cp=sqrt(2)*CycPilot(31);        % get one period of the cyclic preamble; length = 31+1 = 32 
preamble=[cp; cp; cp; cp];  % constrauct the preamble
%
% source
%
info_bits=file2bin('testfile.gif')';      % convert source to bits
%
% generate transmit signal
%
Tb=1e-4 ;               % (nominal) symbol rate
fb=1/Tb;    
epsilon=0;           % fractional error in symbol rate
L=100;                  % oversampling factor
Ts=Tb/L;                % sampling period
fs=1/Ts;
K=8;                    % length of the transmit filter in unit of symbol intervals
fc=1e5;                 % (nominal) carrier frequency
Dfc = 0;
phic = 0;
alpha=0.5;              % transmit filter roll-off factor 

c=1;                    % multipath channel
% c=[1 zeros(1,91) 0.4];
% c=[1 zeros(1,60) 2 zeros(1,123) 0.25]; 
% c=[1 zeros(1,67) 0.75 zeros(1,145) 0.4]; 
%c=[1 zeros(1,75) 0.6 zeros(1,103) 0.2]; 

s=bits2QPSK([1 0 0 0 0 1 0 1 0 1 1 1]);          % convert bits to symbols
% s=[preamble; s];                 % add preamble
pT=sr_cos_p(L*K,L,alpha);        % transmitter filter
xBB=conv(pT,expander(s,L));      % baseband pulse-shaping 

%xBB=conv(pT,expander(s,L));      % baseband pulse-shaping 
% xBB=resample(xBB,epsilon);       % change the sampling rate fs to (1+epsilon)fs
t=(0:length(xBB)-1)'*Ts;         % time
xRF=real(xBB.*exp(1j*(2*pi*(fc+Dfc)*t+phic)));       % modulation
xRF=conv(c,xRF);

yBB = 2*exp(-1j*2*pi*fc*t).*xRF; % desired baseband signal

%filter out out-of-spectrum components
pR=pT; % receiver filter
y = conv(pR, yBB);

figure(2),plot(y,'b')
%plot(y(1:L:end),'b.')
axis('square')
xlabel('real part')
ylabel('imaginary part')
axis([-2 2 -2 2])

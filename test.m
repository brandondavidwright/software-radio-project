%
% QAM Transceiver: 
% This program provide an skeleton for simulation of a quadratue
% amplitude modulated (QAM) transmission system.
%
%%%%%%%%%%%%%%%%%%%%%%
%     PARAMETERS     %
%%%%%%%%%%%%%%%%%%%%%%
Tb=0.0001;      % Symbol/baud period 
L=100;          % Number of samples per symbol period
Ts=Tb/L;        % Sampling period
fc=100000;      % Carrier frequency at the transmitter
delta_c=5;      % Carrier frequency offset
phi_c=0;        % Carrier phase offset
alpha=0.5;      % Roll-off factor for the (square-root) raised cosine filters 
K=8;            % Half of the length of the square-root raised cosine filters in Tb. 
sigma_v=0;      % Standard deviation of channel noise
c=1;            % Channel impulse response

%%%%%%%%%%%%%%%%%%%%%%
%       SOURCE       %
%%%%%%%%%%%%%%%%%%%%%%
N=1000;
b=sign(randn(N,1));

%%%%%%%%%%%%%%%%%%%%%%
%       CODER        %
%%%%%%%%%%%%%%%%%%%%%%
s=b(1:2:end)+i*b(2:2:end);    % 4QAM/QPSK

%%%%%%%%%%%%%%%%%%%%%%
% TRANSMIT FILTERING %
%%%%%%%%%%%%%%%%%%%%%%
pT=sr_cos_p(K,L,alpha);       % Transmit filter
xbbT=conv(expander(s,L),pT);  % Pulse-shaped transmit signal

%%%%%%%%%%%%%%%%%%%%%%
%     MODULATION     %
%%%%%%%%%%%%%%%%%%%%%%
t=[0:length(xbbT)-1]'*Ts;      % Set the time indices
xT=real(exp(i*2*pi*fc*t).*xbbT);

%%%%%%%%%%%%%%%%%%%%%%
%      CHANNEL       %
%%%%%%%%%%%%%%%%%%%%%%
xR=conv(c,xT);
xR=xR+sigma_v*randn(size(xR)); % Received signal

%%%%%%%%%%%%%%%%%%%%%%
%    DEMODULATION    %
%%%%%%%%%%%%%%%%%%%%%%
t=[0:length(xR)-1]'*Ts;         % Set the time indices
xbbR=2*exp(-i*(2*pi*(fc+delta_c)*t-phi_c)).*xR;

%%%%%%%%%%%%%%%%%%%%%%
% RECEIVE FILTERING  %
%%%%%%%%%%%%%%%%%%%%%%
pR=pT;    
y=conv(xbbR,pR);

%%%%%%%%%%%%%%%%%%%%%%
%      DISPLAY       %
%%%%%%%%%%%%%%%%%%%%%%
axes('position',[0.55 0.25 0.4 0.4])
plot(y(15:L:end),'.b')
axis('square')
xlabel('real part')
ylabel('imaginary part')
text(-0.15,-2.2,'(b)')




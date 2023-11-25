%
% Baseband PAM Transceiver: 
% This program provide an skeleton for simulation of a baseband pulse
% amplitude modulated (PAM) transmission system.
%
%%%%%%%%%%%%%%%%%%%%%%
%     PARAMETERS     %
%%%%%%%%%%%%%%%%%%%%%%
Tb=0.0001;      % Symbol/baud period 
L=100;          % Number of samples per symbol period
Ts=Tb/L;        % Sampling period
alpha=0.5;      % Roll-off factor for the (square-root) raised cosine filters 
N=8*L;          % N+1 is the length of the square-root raised-cosine filter. 
sigma_v=0;      % Standard deviation of channel noise
c=1;            % Channel impulse response
%c=[1; zeros(80,1); 0.5; zeros(119,1)];
%%%%%%%%%%%%%%%%%%%%%%
%       SOURCE       %
%%%%%%%%%%%%%%%%%%%%%%
N=1000;
b=sign(randn(N,1));

%%%%%%%%%%%%%%%%%%%%%%
%       CODER        %
%%%%%%%%%%%%%%%%%%%%%%
s=b;                        % 2-level PAM
%s=b(1:2:end)+2*b(2:2:end); % 4-level PAM

%%%%%%%%%%%%%%%%%%%%%%
% TRANSMIT FILTERING %
%%%%%%%%%%%%%%%%%%%%%%
pT=sr_cos_p(N,L,alpha);     % Transmit filter: 
xT=conv(expander(s,L),pT);  % Transmit signal

%%%%%%%%%%%%%%%%%%%%%%
%      CHANNEL       %
%%%%%%%%%%%%%%%%%%%%%%
xR=conv(c,xT);
xR=xR+sigma_v*randn(size(xR)); % Received signal


%%%%%%%%%%%%%%%%%%%%%%
% RECEIVE FILTERING  %
%%%%%%%%%%%%%%%%%%%%%%
pR=pT;    
y=conv(xR,pR);

%%%%%%%%%%%%%%%%%%%%%%
%      DISPLAY       %
%%%%%%%%%%%%%%%%%%%%%%
%
% Here, the eye-diagram of the received PAM signal is displayed.
%
y=reshape(y,2*L,length(y)/(2*L));   % Reshape the received signal 
                                    % vector to columns of duration 2T_b.
y=y(:,8:end-8);     % Delete the transient parts. 
t=[0:1/L:2-1/L];    % Time is unit of symbol period.
axes('position',[0.25 0.25 0.5 0.5])
plot(t,y,'b')   
xlabel('t/T_b')


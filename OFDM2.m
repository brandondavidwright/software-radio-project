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
Los=4;
Ts=Tb/L;      % Sampling period
fc=50000;       % Carrier frequency at the transmitter
Dfc=0;          % Carrier frequency offset
phic=0;         % Carrier phase offset
K=6;            % Length of transmit and receive filters in Tb. 
sigmav=0.000;       % Standard deviation of channel noise
c=1;   % Channel impulse response
%c=[1; zeros(140,1); 0.5; zeros(188,1); 0.3];
%%%%%%%%%%
% Source %
%%%%%%%%%%
b=sign(randn(2*Nactive*Nsymbols,1));
%%%%%%%%%
% coder %
%%%%%%%%%
sfreq=b(1:2:end)+i*b(2:2:end);   % convert bits to QPSK symbols
                                 % These are in the frequency domain.
stime=OFDMmod(sfreq,N,Nactive,Ncp,Los);   % This functions forms the OFDM symbols,
                                    % performs subcarriers modulation (IFFT),
                                    % and add cyclic prefix.
%%%%%%%%%%%%%%%%%%%%%%
% TRANSMIT FILTERING %
%%%%%%%%%%%%%%%%%%%%%%
pT=firpm(30,2*[0 Nactive/2 Los*64-Nactive/2 32*L]/64/L,[1 1 0 0],[1 10])'; 
% alpha=0.25; gamma=1;
% pR=sr_Nyquist_p(K*L,L,alpha,gamma);
xbbT=conv(expander(stime,L/Los),pT);
spec_analysis(xbbT,1)
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
spec_analysis(xR,1)
%%%%%%%%%%%%%%%%%%%%%%
%    DEMODULATION    %
%%%%%%%%%%%%%%%%%%%%%%
t=[0:length(xR)-1]'*Ts;       % Set the time indices
xbbR=2*exp(-i*(2*pi*(fc+Dfc)*t-phic)).*xR;
spec_analysis(xbbR,1)
%%%%%%%%%%%%%%%%%%%%%%
% RECEIVE FILTERING  %
%%%%%%%%%%%%%%%%%%%%%%
pR=pT;  
% alpha=0.25; gamma=1;
% pR=sr_Nyquist_p(K*L,L,alpha,gamma);
y=conv(xbbR,pR);
spec_analysis(y,1)
y=y(1:L:end);                   % Decimation


h=zeros(32,1);
x=stime(1:Los:end);
mu=0.1;
for n=32:length(x)
    xtdl=x(n:-1:n-31);
    e=(y(n-3)-h'*xtdl)';
    h=h+mu/(xtdl'*xtdl)*xtdl*e;
    xi(n)=abs(e)^2;
end



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Remove CP, apply DFT and FEQ %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% n=Ncp;
% w=ones(Nactive,1); % mu=0.1;
% m=1; mu=0.1; shat=[];
% p=conv(pT,pR); c=c.*exp(-j*2*pi*[0:length(c)-1]'*Ts*fc); cBB=conv(c,p);
% cBB=cBB(1:L:end); %cBB=[cBB(2:end);cBB(1)]; 
% CBB=fft(cBB,N); 
% k=1;
% while n<=length(y)-N
%     A1=fft(y(n+1:n+N));
%     %A1=A1./CBB;         % Equalization
%     A1=[A1(2:Nactive/2+1); A1(end-Nactive/2+1:end)];
%     e=conj(sfreq(m:m+Nactive-1)-conj(w).*A1);
%     w=w+2*mu*A1.*(abs(A1).^(-2)).*e;
%     %x=[x; conj(w).*A1];
%     shat=[shat; A1];
%     n=n+Ncp+N;
%     m=m+Nactive;
%     figure(3),plot(conj(w).*A1,'.'),axis([-1.5 1.5 -1.5 1.5]),pause(0.1)
%     xi(k)=real(e'*e)/52; k=k+1;
% %    figure(2),plot(A1,'.'),axis([-2.5 2.5 -2.5 2.5]),pause(0.1)
% end
figure,semilogy(xi)

figure,plot(abs(h))

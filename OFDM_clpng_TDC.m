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
L=40;        %%%%% Interpolation factor
M=L;        %%%%% Add raised cosine roll-off.
stime=OFDMmod_wndng(sfreq,N,Nactive,Ncp,M,L);   % This functions forms the OFDM symbols,
                                    % performs subcarriers modulation (IFFT),
                                    % and add cyclic prefix.
stime=stime;
%%%%%%%%%%%%%%%%%%%%%%
% TRANSMIT FILTERING %
%%%%%%%%%%%%%%%%%%%%%%
% pT=firpm(K*L,[0 24/32/L 40/32/L 1],[1 1 0 0],[1 10])'; 
 
xbbT=stime;
stdxbbT=std(xbbT);
cmp=abs(xbbT)<(2*stdxbbT);
xclipped=xbbT.*cmp+2*stdxbbT*(1-cmp).*xbbT./(abs(xbbT)+1e-6);
xd=abs(xbbT-xclipped);
List=[];i1=0;i2=0;
for n=2:length(xd)
    if (xd(n-1)<1e-10)&(xd(n)>1e-10)
        i1=n-1;
    elseif (xd(n-1)>1e-10)&(xd(n)<1e-10)
        i2=n;
    end
    if (i1~=0)&(i2~=0)
        List=[List round((i1+i2)/2)];i1=0;i2=0;
    end
end

g=prolate(2*L,0.05); % prolate pulse for clipping compensation.
d=zeros(size(xd));
for n=List
    d(n-L:n+L)=(abs(xbbT(n))-2*stdxbbT)*g.*xbbT(n-L:n+L)./abs(xbbT(n-L:n+L));
end
xmclipped=xbbT-d;
xclipped1=xclipped;

cmp=abs(xmclipped)<(2*stdxbbT);
xclipped=xmclipped.*cmp+2*stdxbbT*(1-cmp).*xmclipped./(abs(xmclipped)+1e-6);
xd=abs(xmclipped-xclipped);
List=[];i1=0;i2=0;
for n=2:length(xd)
    if (xd(n-1)<1e-10)&(xd(n)>1e-10)
        i1=n-1;
    elseif (xd(n-1)>1e-10)&(xd(n)<1e-10)
        i2=n;
    end
    if (i1~=0)&(i2~=0)
        List=[List round((i1+i2)/2)];i1=0;i2=0;
    end
end

d=zeros(size(xd));
for n=List
    d(n-L:n+L)=(abs(xmclipped(n))-2*stdxbbT)*g.*xmclipped(n-L:n+L)./abs(xmclipped(n-L:n+L));
end
xmclipped=xmclipped-d;

cmp=abs(xmclipped)<(2*stdxbbT);
xmclipped=xmclipped.*cmp+2*stdxbbT*(1-cmp).*xmclipped./(abs(xmclipped)+1e-6);

xclipped=xclipped1;



XX=0;XXclipped=0;XXmclipped=0;
K=0;
N=(N+Ncp)*L+M;
for k=1:N:length(xbbT)-N
    [X,F]=spec_analysis([xbbT(k:k+N-1)],1);
    XX=XX+X.^2;
    [X,F]=spec_analysis([xclipped(k:k+N-1)],1);
    XXclipped=XXclipped+X.^2;
    [X,F]=spec_analysis([xmclipped(k:k+N-1)],1);
    XXmclipped=XXmclipped+X.^2;
    K=K+1;
end
XXclipped=10*log10(XXclipped/max(XXclipped));
XXmclipped=10*log10(XXmclipped/max(XXmclipped));
XX=10*log10(XX/max(XX));
figure,axes('position',[0.25 0.25 0.5 0.5])
plot(F,XX,'k',F,XXclipped,'r',F,XXmclipped,'g')
xlabel('Normalized Frequency, $f$','interpreter','latex')
ylabel('Amplitude, dB')
  
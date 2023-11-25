%
% TR_ELG: Timing recovery using early-late gating.
%
Tb=0.0001; L=100; M1=20; Ts=Tb/L; fs=1/Ts; fc=100000; 
Dfc=10; N=8*L; phic=pi/4; sigmav=0; alpha=0.5; c=1; 
%c = [0.7; zeros(93,1);  1; zeros(107,1); -0.4];
Nb=1200;b=sign(randn(Nb,1)); 
M=input('QAM size (4, 16, 64, 256) =');
if M==4 s=b(1:2:end)+i*b(2:2:end);    
elseif M==16 s=2*b(1:4:end)+b(2:4:end)+i*(2*b(3:4:end)+b(4:4:end));
elseif M==64 s=4*b(1:6:end)+2*b(2:6:end)+b(3:6:end)+...
        j*(4*b(4:6:end)+2*b(5:6:end)+b(6:6:end));
elseif M==256 s=8*b(1:8:end)+4*b(2:8:end)+2*b(3:8:end)+b(4:8:end)+...
        j*(8*b(5:8:end)+4*b(6:8:end)+2*b(7:8:end)+b(8:8:end));   
else print('Error! M should be 4, 16, 64 or 256'); end 
pT=sr_cos_p(N,L,alpha); xbbT=conv(expander(s,L),pT);  
t=[0:length(xbbT)-1]'*Ts; xT=real(exp(i*2*pi*fc*t).*xbbT);
xR=conv(c,xT); xR=xR+sigmav*randn(size(xR)); 
t=[0:length(xR)-1]'*Ts; y=2*exp(-i*(2*pi*(fc-Dfc)*t-phic)).*xR;
pR=pT; y=conv(y,pR); y=y(1:M1:end); fs1=fs/M1;y4=y.^4;
[X,F]=spec_analysis(y4,fs1);
figure, axes('position',[0.1 0.25 0.8 0.5]), plot(F,X,'k')
xlabel('FREQUENCY, Hz'), ylabel('AMPLITUDE')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                    %
%           EXTENSITONS              %
%                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extention 1: Coarse carrier acquisition    %
% Comment out this extension, if you wish    %
% to look at the fine acquisition algorithm. %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [xmax,imax]=max(X);  Dfc_est=F(imax)/4;
% y1=y.*exp(-j*2*pi*Dfc_est*Ts*M1*[0:length(y)-1]');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Fine carrier acquisition and tracking    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Extension 2: lowpass filtering and decimation
%
Dfcmax=100; M2=25; fs2=fs1/M2; fpb=4*Dfcmax; fsb=fs2-fpb; 
Ng1=124; g1=firpm(Ng1,[0 fpb fsb fs1/2]/(fs1/2),[1 1 0 0]);
beta=filter(g1,1,y4);  beta=beta(1:M2:end); 
[X,F]=spec_analysis(beta,fs2);
figure, axes('position',[0.1 0.25 0.8 0.5]), plot(F,X,'k')
xlabel('FREQUENCY, Hz'), ylabel('AMPLITUDE')
%
% Extension 3: adaptive line enhancer (bandpass filtering)
%
Nf=100;mu=0.01;f=zeros(Nf,1); eta=zeros(size(beta));
for n=Nf+1:length(beta)
    xtdl=beta(n-1:-1:n-Nf);
    eta(n)=f'*xtdl;
    e=beta(n)-eta(n);
    f=f+2*mu*e'*xtdl/(xtdl'*xtdl);
end
[X,F]=spec_analysis(eta,fs2);
figure, axes('position',[0.1 0.25 0.8 0.5]), plot(F,X,'k')
xlabel('FREQUENCY, Hz'), ylabel('AMPLITUDE')
%
% Extension 4: PLL
%
Ts=1/fs2; mu=Ts;
BL=10;zeta=0.707;wn=2*BL/(zeta+1/(4*zeta));
KL=wn^2*Ts+2*zeta*wn;
alpha=-2*zeta/(2*zeta+wn*Ts);
c1=0;e1=0;phi=zeros(size(eta));eta=-eta;
for n=2:length(eta)
    e=phasedetector(eta(n-1),exp(j*phi(n-1)*4),e1);
    c=(KL/4)*mu*(e+alpha*e1)+c1;
    phi(n)=c+phi(n-1);
    if abs(phi(n))>pi
        phi(n)=phi(n)-2*pi*sign(phi(n));
    end
    c1=c;e1=e;    % values for the next iteration
end
%
% Extension 5: Interpolation and frequency/phase compensation.
%
theta=expander(phi,M2);                                % Expansion
for n=1:length(phi)-1
    a=theta(M2*n+1); b=theta(M2*(n-1)+1);
    amb=a-b;
    if abs(amb+2*pi)<abs(amb)
        a=a+2*pi;
    elseif abs(amb-2*pi)<abs(amb)
        a=a-2*pi;
    else
    end
    theta(M2*(n-1)+2:M2*n)=b+([1:M2-1]'/M2)*(a-b);
end
Delta=Ng1/2; y=[zeros(Delta,1); y(1:length(theta)-Delta)];
y1=y.*exp(-j*(theta));          % Carrier and phase compensated output.
%
% Generate a movie of eye-patterns 
%
figure
k=3+M2;                          % Timing phase
while (k+1000)<length(y1)    % Two-dimensional eye-pattern of the phase corrected 
    plot(y1(k-M2:5:k+1000),'.'),hold on,plot([-1.5 1.5],[-1.5,1.5]),pause(1),hold off   % received signal.
    k=k+1000;
end

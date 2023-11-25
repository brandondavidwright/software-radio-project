%
% CRExp1.m
%
Tb=0.0001; L=100; Ts=Tb/L; fc=100000; fs=1/Ts;  
alpha=0.5;  N=8*L;  sigma_v=0;  c=1;      
s=sign(randn(1000,1));
pT=sr_cos_p(N,L,alpha); xT=conv(expander(s,L),pT);  
t=[0:length(xT)-1]'*Ts; xT=cos(2*pi*fc*t).*xT;
xR=conv(c,xT); xR=xR+sigma_v*randn(size(xR));
xR2=xR.^2;
[X,F]=spec_analysis(xR2,fs);
figure(2),axes('position',[0.1 0.25 0.8 0.5])
plot(F,X,'k')
xlabel('FREQUENCY, Hz')
ylabel('AMPLITUDE')

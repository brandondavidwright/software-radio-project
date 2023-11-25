%
% CREp2.
%
Tb=0.0001; L=100; M1=20; Ts=Tb/L; fs=1/Ts; fc=100000; 
Dfc=10; N=8*L; phic=0.5; sigma_v=0; alpha=0.5; c=1; 
c = [0.7; zeros(93,1);  1; zeros(107,1); -0.4];
b=sign(randn(12000,1)); 
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
xR=conv(c,xT); xR=xR+sigma_v*randn(size(xR)); 
t=[0:length(xR)-1]'*Ts; y=2*exp(-i*(2*pi*(fc-Dfc)*t-phic)).*xR;
pR=pT; y=conv(y,pR); y=y(1:M1:end); fs1=fs/M1;
[X,F]=spec_analysis(y.^4,fs1);
figure, axes('position',[0.1 0.25 0.8 0.5]), plot(F,X,'k')
xlabel('FREQUENCY, Hz'), ylabel('AMPLITUDE')

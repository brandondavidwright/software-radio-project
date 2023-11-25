%
% DDCR.m
%
Tb=0.0001; L=100; M1=20; Ts=Tb/L; fs=1/Ts; fc=100000; 
Dfc=0; N=8*L; phic=pi/8; sigmav=0.05; alpha=0.5; c=1; 
b=sign(randn(1000,1)); 
M=4; s=b(1:2:end)+i*b(2:2:end);    
pT=sr_cos_p(N,L,alpha); xbbT=conv(expander(s,L),pT);  
t=[0:length(xbbT)-1]'*Ts; xT=real(exp(i*2*pi*fc*t).*xbbT);
xR=conv(c,xT); xR=xR+sigmav*randn(size(xR)); 
t=[0:length(xR)-1]'*Ts; y=2*exp(-i*(2*pi*(fc-Dfc)*t-phic)).*xR;
pR=pT; y=conv(y,pR); y=y(1:L:end);y=y(9:end-8); 
phi=zeros(size(y)); s1=zeros(size(y)); mu=0.01;
for n=1:length(y)-1
    s1(n)=y(n)*exp(-j*phi(n));
    s2=sign(real(s1(n)))+j*sign(imag(s1(n)));
    s12=s1(n)*s2'; e=imag(s12)/real(s12);
    phi(n+1)=phi(n)+mu*e;
end
figure(1),axes('position',[0.25 0.25 0.5 0.5]),plot(phi,'k')
xlabel('n'),ylabel('\phi[n]')
figure(2),axes('position',[0.25 0.25 0.5 0.5]),
plot(s1(1:end-1),'*k'),axis('square'),hold on,
plot([-2 2],[-2 2],'k--',[-2 2],[2 -2],'k--'),
xlabel('Real part'),ylabel('Imaginary part')
hold off
% plot(y,'ok')
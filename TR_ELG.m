%
% TR_ELG: Timing recovery using early-late gating.
%
Tb=0.0001; L=100; Ts=Tb/L; fs=1/Ts; fc=100000; 
Dfc=0; N=8*L; phic=0; sigma_v=0; alpha=0.5; c=1; 
b=sign(randn(4000,1)); 
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
t=[0:length(xR)-1]'*Ts; y=exp(-i*(2*pi*(fc-Dfc)*t-phic)).*xR;
pR=pT; y=conv(y,pR); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  TIMING RECOVER: Early-late Gating    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta=0;               
mu0=0.01;dtau=12;mu=mu0*(L/4)/dtau;Ly=length(y);
kk=1;yp=0;ym=0;start=5*L+1;tau=0.3*ones(1,floor((Ly-start)/L));
for k=start:L:length(tau)*L
    tauTb=round(tau(kk)*L);
    yp=sqrt(1-beta^2)*y(k+tauTb+dtau)-beta*yp;
    ym=sqrt(1-beta^2)*y(k+tauTb-dtau)-beta*ym;
    tau(kk+1)=tau(kk)+mu*(abs(yp)^2-abs(ym)^2);
    kk=kk+1;
end
figure, axes('position',[0.1 0.25 0.8 0.5]), plot(tau(1:kk),'k')
xlabel('Iteration Number, n'), ylabel('\tau[n]')

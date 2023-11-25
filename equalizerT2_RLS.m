%
% Tb/2-spaced equalizer (RLS algorithm)
%
clear all, close all
Tb=0.0001; L=100; Ts=Tb/L; fs=1/Ts; fc=100000; 
Dfc=0; N=10*L; phic=0; sigma_v=0.01/2; alpha=0.5; 
TmgPhase=51; lambda=0.99;
c=[1 zeros(1,91) 0.4]';
% c=[0.5 zeros(1,60) 1 zeros(1,123) 0.25]'; 
% c=[1 zeros(1,67) 0.75 zeros(1,145) 0.4]'; 
% c=[1 zeros(1,75) 0.6 zeros(1,103) 0.2]';
c=c/sqrt(c'*c);
b=sign(randn(2000,1)); 
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
pR=pT; y=conv(y,pR); 
y=y(TmgPhase:L/2:end);
xi=zeros(length(y)/2,1); x=xi;
N=31;           % Equalizer order
N1=N+1;          % Equalizer length
Delta=round((N1/2+(length(c)+length(pT))/L)/2);
w=zeros(N1,1); nn=round((N+Delta)/2);nn0=nn
Psi_inv=10000*eye(N1); 
for n=N1+Delta:2:length(y)-Delta-N1
    tdl=y(n:-1:n-N1+1);
    u=Psi_inv*tdl;
    k=u/(lambda+tdl'*u);
    x(nn)=w'*tdl;
    e=s(nn-Delta)-x(nn);
    w=w+k*e';
    Psi_inv=(1/lambda)*(Psi_inv-k*(tdl'*Psi_inv));
    xi(nn-nn0+1)=abs(e)^2;
    if rem(nn,400)==0 
        figure(1),plot(x(nn-398:nn),'.'),pause(0.1) 
    end
    nn=nn+1;
end
figure(2),semilogy(xi)

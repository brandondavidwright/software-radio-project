%
% MATLAB Script MLPLLtest.m
%
% Main
%
fc=100; Ts=0.001; sigmav=sqrt(0.01); a=1; KL=10;
t=[0:2000]*Ts; theta=40*pi/41*ones(size(t));
phi1=MLPLL1(fc,Ts,a,KL,t,theta,sigmav);
phi1a=MLPLL1a(fc,Ts,a,KL,t,theta,sigmav);
phimean=MLPLLmean(fc,Ts,KL,a,t,theta);
axes('position',[0.25 0.25 0.5 0.5])
plot(t,phi1,'--k',t,phi1a,':k',t,phimean,'-k',t,theta,'LineWidth',1);
xlabel('time, sec.')
ylabel('\phi , radian')
legend('(8.88)','(8.89)','(8.82)')

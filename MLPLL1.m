%
% ML PLL based on (\ref{eqn:ch8:LMS2a})\\
%
function phi=MLPLL1(fc,Ts,a,KL,t,theta,sigmav)
x=a*cos(2*pi*fc*t+theta)+sigmav*randn(size(t));
phi=zeros(size(t)); y=phi; mu=Ts;
for n=1:length(t)-1
    y(n)=a*cos(2*pi*(n-1)*fc*Ts+phi(n));
    phi(n+1)=phi(n)-2*(mu*KL/a)*sin(2*pi*(n-1)*fc*Ts+phi(n))*(x(n)-y(n));
end
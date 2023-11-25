%
% Mean PLL based on  (\ref{eqn:ch8:LMS5})
%
function phi=MLPLLmean(fc,Ts,KL,a,t,theta)
phi=zeros(size(t)); mu=Ts;
for n=1:length(t)-1
    phi(n+1)=phi(n)+(mu*KL)*sin(theta(n)-phi(n));
end

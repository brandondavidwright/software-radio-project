clear all
n_ensemble=10000; Ldata=300; sigman=sqrt(0.001);
h=[0.2 1 -0.1]'+i*[-0.3 0.5 -0.7]'; N=length(h)-1;
mu=0.01;
xi=zeros(Ldata,1);
for k=1:n_ensemble
    w=zeros(size(h));
    x=sign(randn(Ldata,1))+i*sign(randn(Ldata,1));
    d=conv(h,x);
    d=d(3:end)+(sigman/sqrt(2))*(randn(Ldata,1)+i*randn(Ldata,1));
    for n=N+1:Ldata;
        xtdl=x(n:-1:n-N);
        e=d(n-N)-w'*xtdl;
        w=w+2*mu*e'*xtdl;
        xi(n-N)=xi(n-N)+abs(e)^2;
    end
end
xi=xi/n_ensemble;
n=[0:Ldata-1]';
axes('position',[0.25 0.25 0.5 0.5])
semilogy(n,xi,'-k')
axis([0 Ldata 0.001 10])
xlabel('Number of iterations')
ylabel('\xi')

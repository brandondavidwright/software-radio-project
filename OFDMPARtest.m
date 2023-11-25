Np=1000;
PAR=zeros(Np,1);
Na=26;N=64;
for k=1:Np
    X=[0; sign(randn(Na,1))+j*sign(randn(Na,1)); zeros(N-2*Na-1,1); sign(randn(Na,1))+j*sign(randn(Na,1))];
    % x=length(X)*ifft(X);
    % figure(1), plot(abs(x))
    Y=[X(1:Na+1); zeros(N*7,1); X(Na+2:end)];
    plot(Y,'.'),hold on
    y=length(Y)*ifft(Y);
    
    % figure(2), plot(abs(y)),hold
    stdy=std(y);
    cmp=abs(y)<(2.5*stdy);
    y=y.*cmp+2.5*stdy*(1-cmp).*y./(abs(y)+1e-6);
    % plot(abs(y)./stdy),pause
%     for n=1:length(y)
%         if abs(y(n))>3*std(y)
%             y(n)=3*stdy*y(n)/abs(y(n));
%         end
%     end
    Y=fft(y).*[ones(N/2,1); zeros(7*N,1); ones(N/2,1)];
    plot(Y/500,'.r'), hold off, pause(0.1)
    y=ifft(Y); % plot(abs(y),'r'),pause
    % figure(3), plot(abs(fft(y,1024)))
    % [max(abs(x)) max(abs(y)) 20*log10(max(abs(y))/std(y))]
    PAR(k)=20*log10(max(abs(y))/stdy);
end
a=sort(PAR);
semilogy(a,[Np:-1:1]/Np)

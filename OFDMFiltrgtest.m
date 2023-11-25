
Na=26*4;N=256;
X=[0; sign(randn(Na,1))+j*sign(randn(Na,1)); zeros(N-2*Na-1,1); sign(randn(Na,1))+j*sign(randn(Na,1))];
x=length(X)*ifft(X);
figure(1), plot([0:length(x)-1],abs(x)),hold
Y=[X(1:Na+1); zeros(N*7,1); X(Na+2:end)];
y=length(Y)*ifft(Y); 
figure(1), plot([0:length(y)-1]/8,abs(y),'r')
figure(3), plot(abs(fft(y,1024*4)))
% [max(abs(x)) max(abs(y)) 20*log10(max(abs(y))/std(y))]

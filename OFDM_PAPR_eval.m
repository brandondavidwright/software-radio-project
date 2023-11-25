Np=100000; PAPR=zeros(Np,1);
Nactive=26;N=64;
for k=1:Np
            % Generate a set of QPSK symbols for one OFDM symbol
    X=[0; sign(randn(Na,1))+j*sign(randn(Na,1)); ...
        zeros(N-2*Na-1,1); sign(randn(Na,1))+j*sign(randn(Na,1))];
    Y=[X(1:Na+1); zeros(N*7,1); X(Na+2:end)]; % Interpolate and
    y=ifft(Y);                          % convert to time domain
    PAPR(k)=20*log10(max(abs(y))/std(y));
end
figure,axes('position',[0.25 0.25 0.5 0.5])
semilogy(sort(PAPR),[Np:-1:1]/Np)   % Plot CCDF
axis([0 12 1e-4 1])
xlabel('PAPR, dB'),ylabel('CCDF') 
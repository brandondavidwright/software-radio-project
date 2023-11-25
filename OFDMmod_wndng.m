%
% Modulator part of an OFDM transceiver
%
function stime=OFDMmod(sfreq,N,Nactive,Ncp,M,L)
stime=[];
for k=1:Nactive:length(sfreq)
    A1=sfreq(k:k+Nactive-1);    % Take a block of QAM symbols
    A=zeros(N*L,1);               % and put them at right places in 
    A(2:Nactive/2+1)=A1(1:Nactive/2);   % the vector A, before applying
    A(end-Nactive/2+1:end)=A1(Nactive/2+1:end); % IFFT.
    a=ifft(A);              % modulate/convert to time-domain
    g=[(1-cos(pi*([-L*Ncp:-L*Ncp+M-1]+L*Ncp+1)/M))/2 ones(1,L*N-M+L*Ncp) (1+cos(pi*([L*N-M:L*N-1]-L*N+M)/M))/2]';
    a=[a((N-Ncp)*L+1-M:end); a].*g;    % add CP
    stime=[stime; a];       % concatenate symbols
end

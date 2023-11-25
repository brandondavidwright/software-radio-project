%
% This program evaluates the minimum mean-squared error (MMSE) of the symbol-spaced 
% and half symbol-spaced equalizers. The equivalent baseband channel cBB(t)
% is obtained according to (3.54). The MMSE values are then calculated, following 
% the formulations in Section 11.3.1. 
%
% The variance of data symbols s[n], sigmas2, is set equal to one.
% The noise variance sigmanuc2 is also set equal to one.  
%
clear all
Tb=0.0001; L=100; Ts=Tb/L; fs=1/Ts; fc=100000; 
alpha=0.5; sigmas=1; sigmanuc=0.01; 
c=[1 zeros(1,91) 0.4];
c=[0.5 zeros(1,60) 1 zeros(1,123) 0.25]; 
c=[1 zeros(1,67) 0.75 zeros(1,145) 0.4]; 
c=[1 zeros(1,75) 0.6 zeros(1,103) 0.2]; 
c=c/sqrt(c*c');
pT=sr_cos_p(16*L,L,alpha)'; 
pR=pT; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Construction of the equivalent baseband channel    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p=conv(pT,pR);
c=c.*exp(-j*2*pi*[0:length(c)-1]*Ts*fc);
cBB=conv(c,p);
pR=sqrt(L/2)*pR(1:L/2:end);

N=31;           % Equalizer order
N1=N+1;          % Equalizer length

%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Tb spaced equalizer   %
%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=1:L
    cBB0=cBB(k:L:end);
    Delta=round((length(cBB(1:L:end))+N1)/2);
    C=toeplitz([cBB0 zeros(1,N)],[cBB0(1) zeros(1,N)]);
    P0=toeplitz([pR(1:2:end) zeros(1,N)],[pR(1) zeros(1,N)]);
    P1=toeplitz([pR(2:2:end) zeros(1,N)],[pR(2) zeros(1,N)]);
    Q=[C; (sigmanuc/sigmas)*P0; (sigmanuc/sigmas)*P1];
    d=[zeros(Delta,1); 1;zeros(length(Q(:,1))-(Delta+1),1)];
    Ryy=(Q'*Q+1e-12*eye(N1));
    pyd=(Q'*d);
    w=Ryy\pyd;
    mmse0(k)=sigmas^2*real(1-w'*pyd);
    spower(k)=sum(abs(cBB0).^2);
    D(k)=Delta;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Tb/2 spaced equalizer  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
cBBf=cBB(1:L/2:end);Delta=round((length(cBBf)+N1)/4); 
for k=1:L
    cBBf=cBB(k:L/2:end);
    C=toeplitz([cBBf zeros(1,N)],[cBBf(1) zeros(1,N)]); C=C(1:2:end,:);
    P=toeplitz([pR zeros(1,N)],[pR(1) zeros(1,N)]);
    Q=[C; (sigmanuc/sigmas)*P];
    d=[zeros(Delta,1); 1;zeros(length(Q(:,1))-(Delta+1),1)];
    Ryy=(Q'*Q+1e-6*eye(N1));
    pyd=(Q'*d);
    w=Ryy\pyd;
    mmsef1(k)=sigmas^2*real(1-w'*pyd);
end

N=61;
N1=N+1;
cBBf=cBB(1:L/2:end);Delta=round((length(cBBf)+N1)/4);
for k=1:L
    cBBf=cBB(k:L/2:end);
    C=toeplitz([cBBf zeros(1,N)],[cBBf(1) zeros(1,N)]); C=C(1:2:end,:);
    P=toeplitz([pR zeros(1,N)],[pR(1) zeros(1,N)]);
    Q=[C; (sigmanuc/sigmas)*P];
    d=[zeros(Delta,1); 1;zeros(length(Q(:,1))-(Delta+1),1)];
    Ryy=(Q'*Q+1e-6*eye(N1));
    pyd=(Q'*d);
    w=Ryy\pyd;
    mmsef2(k)=sigmas^2*real(1-w'*pyd);
end

figure
subplot(2,1,1),plot([0:L-1]/L,spower); ylabel('Signal Power')
subplot(2,1,2),semilogy([0:L-1]/L,mmse0,'-',[0:L-1]/L,mmsef1,'--',[0:L-1]/L,mmsef2,'-.')
xlabel('Timing Phase'),ylabel('MMSE'),legend('T_b spaced equalizer (N=31)','T_b/2 spaced equalizer (N=31)','T_b/2 spaced equalizer (N=61)')


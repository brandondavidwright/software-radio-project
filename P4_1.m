%
% P4-1: This program is designed to help the reader to apprecite the impact
% of anti-aliasing filter.
%

%
% The first part of the program samples the signal without using any
% anti-aliasing filter and reconstruct the signal at a much higher rate.
% The recontructed signal is then compared with the continuous time signal
% which is presented by a dense set of samples.
%
clear all
Ts=0.25;L=10;
td=[-10:Ts:10];
xd=exp(-abs(td));
tc=[-10:Ts/L:10];
xc=exp(-abs(tc));
xr=zeros(1,length(xc));
N=(length(xd)-1)/2;
for n=-N:N
    nn=n+N+1;
    xr=xr+xd(nn)*sinc((tc-n*Ts)/Ts);
end
figure(1),plot(tc,xc,tc,xr)
e=xc-xr; e2=e*e'

%
% Here, anti-aliasing filter is applied before sampling. As above
% continuous-time signal is approximated by a dense grid of samples.

h=firpm(100,[0 1/L 1.5/L 1],[1 1 0 0]);
xf=conv(h,xc);
xf=xf(51:end-50);
xfd=xf(1:L:end);
xfr=zeros(1,length(xc));
N=(length(xfd)-1)/2;
for n=-N:N
    nn=n+N+1;
    xfr=xfr+xfd(nn)*sinc((tc-n*Ts)/Ts);
end
figure(2),plot(tc,xc,tc,xfr)
ef=xc-xfr; ef2=ef*ef'
%
% RAISED-COSINE PULSE: h=r_cos_p(N,L,alpha)
% This function generates a raised cosine pulse of length N+1. 
% There are L samples per symbol period.
% alpha is the roll-off factor.
%
function h=r_cos_p(N,L,alpha)
t=[-N/2:1:N/2]/L;
h=zeros(size(t));
for k=1:length(t)
    if (1-(4*alpha^2)*t(k)^2)==0
        h(k)=pi*alpha*sin(pi*alpha*t(k))/((8*alpha^2)*t(k));
    else
        h(k)=cos(pi*alpha*t(k))/(1-(4*alpha^2)*t(k)^2);
    end
end
h=h.*sinc(t);
h=h';
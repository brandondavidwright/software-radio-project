% PLL with extended lock range
% This follows the procedure mentioned in Section 8.4 of the text.
% See Figure 8.14 for some details.
%
% Variables:
%   eta: PLL input 
%   phi: PLL ouput
%   e: phase error
%   K_L and alpha: are the loop filter parameters
%   mu: the integrator parameter
%
fs2=2000;
mu=1/fs2;KL=284.5;alpha=-0.9375;
c1=0;e1=0;eta=exp(j*[0:100]);
phi=zeros(size(eta));
for k=2:length(eta)
    e=phasedetector(eta(k),exp(j*phi(k-1)),e1);
    c=KL*(e+alpha*e1)+c1;
    phi(k)=mu*c+phi(k-1);
    c1=c;e1=e;
end
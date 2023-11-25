%
% Isotropic pulse design
%
% h=isotropic_pf(N,M,alpha);
%
% It designs a square-root raised-cosine pulse-shape with the following
% parameters:
%   N: filter order (filter length = N+1)
%   M: number of samples per symbol period Tb
%   alpha: roll-off factor (between 0 and 1)


function out = isotropic_pf(N, M, alpha)
L = 256; %Sampling rate
K = 6; %Length of the anti-casual filter
t = sqrt(1+alpha)*(-K:1/L:K)';
nIter = 100;

% Number of Hermite Polynomials used for the design
NUM = 4;

% The hermite polynomials
D0 = @(t)exp(-pi*t.^2); 
D4 = @(t)16*exp(-pi*t.^2).*(4*pi^2*t.^4-6*pi*t.^2+3/4);
D8 = @(t)256*exp(-pi*t.^2).*(16*pi^4*t.^8-112*pi^3*t.^6+210*pi^2*t.^4-105*pi*t.^2+105/16);
D12 = @(t)4096*exp(-pi*t.^2).*(64*pi^6*t.^12-1056*pi^5*t.^10+5940*pi^4*t.^8-...
    13860*pi^3*t.^6+51975/4*pi^2*t.^4-31185/8*pi*t.^2+10395/64);
d0 = D0(t);
d4 = D4(t);
d8 = D8(t);
d12 = D12(t);

H = [d0 d4 d8 d12]; 
A= {[0 0], [1 0], [1 1] [2 0]}; %The set of points that we are interested in
lambda = [1 40 40 40]; % Weighting function used
u = [1 0 0 0]'; % desired values at each point

% Computing the cross ambiguity function
G = cell(size(A));
weights = zeros(size(A));

for ii = 1:length(A)
    T = A{ii}(1);
    F = A{ii}(2);
    T_Offset = abs(T*L);
    F_Offset = abs(F*L)*sqrt(1+alpha);
    eexp = exp(-1j*2*pi*t*F_Offset/L);
    H1 = [zeros(T_Offset, NUM); H];
    HS = [d0.*eexp d4.*eexp d8.*eexp d12.*eexp];
    H2 = [HS; zeros(T_Offset, NUM)];
    G{ii} = H2'*H1;
    weights(ii) = lambda(ii);
end

% Initialization
a = [1; zeros(NUM-1,1)];

% Iterative Opitmization
gamma2 = diag(weights.^2);
for i = 1:nIter
    for j = 1:length(G)
        D(j, :) = a'*G{j};
    end
    a_n = (D'*gamma2*D)^-1*D'*gamma2*u;
    a = real(a_n+a)/2;
end

% Producing the final samples
K = (N/2/M);
t = sqrt(1+alpha)*(-K:1/M:K)';
d0 = D0(t);
d4 = D4(t);
d8 = D8(t);
d12 = D12(t);
H = [d0 d4 d8 d12];  
h = H*a;
out = h/sqrt(max(conv(h, h)));

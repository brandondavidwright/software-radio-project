close all; clear;

load('xRF9.mat')

%----------------------------
%---begin part 1-------------
%----------------------------

% examine spectrum of xRF
figure
spec_analysis(xRF, fs)
title("Received signal");

% convert to baseband QAM signal
t = (0:length(xRF)-1)'*Ts; %time of xRF
N = length(cp);

xBB = 2*exp(-1j*2*pi*fc*t).*xRF; % desired baseband signal

% examine spectrum of unfiltered baseband
figure
spec_analysis(xBB,fs);
title("Unfiltered baseband signal");

%filter out out-of-spectrum components
pR=pT; % receiver filter
xBBf = conv(xBB, pR);

figure
spec_analysis(xBBf,fs);
title("filtered baseband signal");
% examine spectrum of unfiltered baseband

% look at eye pattern
figure
eye_pattern(xBBf);
title("eye pattern of cBB")

%----------------------------
%---part 1 to be continued---
%----------------------------
%---begin part 5-------------
%----------------------------

% track symbol rate to recover timing
M = 4;
ytau = ddtr(xBBf, L, M);
xBBd = ytau;

figure
eye_pattern(xBBd)
title("eye pattern of xBBd")

%----------------------------
%---end of part 5------------
%----------------------------
%---part 1 continues---------
%----------------------------

% identify preamble and extract payload
%detect preamble
% extract data symbols from payload and convert to bit stream
preamble = [cp; cp; cp; cp];
payload_start = find_payload_start(xBBd, cp);

payload = xBBd(payload_start:end);
figure
eye_pattern(payload)
title("eye pattern of payload")

bits = QPSK2bits(payload);

% save bits to file
bin2file(bits', "transmitted_file.txt");

%----------------------------
%---end of part 1------------
%----------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Find start of payload of signal y and pilot cp                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function index = find_payload_start(y, cp)
    correlations = zeros(1,length(y));
    for i = 1:1:length(y)
        if i+length(cp)-1 < length(y)
            % find correlation of wind of signal y with pilot cp
            correlations(i) = corr(real(y(i:i+length(cp)-1)), real(cp));
        end
    end
    % Find beginning of last correlated pilot segment
    start_last_cp = find(correlations>0.9, 1, "last");
    % Find beginning of payload
    index = start_last_cp+length(cp);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Track symbol phase of signal y with dec. factor L and order M   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ytau = ddtr(y, L, M)
    mu = 0.05;
    Ly = length(y);
    start = 5*L+1;
    kk = 1;
    tau = 0.3*ones(floor((Ly-start)/L),1);
    dtau = 6; % TODO what is dtau?
    ytau = zeros(size(tau));
    
    % decision directed timing phase recovery method (Section 10.3.2)
    for k = start:L:length(tau)*L-L
        tauTb = round(tau(kk)*L);
        ytau(kk) = y(k + tauTb);
        sk = slicer(y(k+tauTb),M);
        tau(kk+1) = tau(kk) + mu*real((sk-y(k+tauTb)) ...
                *(y(k+tauTb+dtau)-y(k+tauTb-dtau))');
        kk = kk + 1;
    end
end

function eye_pattern(y)
    plot(y, "b")
    xlabel("real part")
    ylabel("imaginary part");
end
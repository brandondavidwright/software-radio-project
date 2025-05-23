close all; clear;

load('xRF1.mat')

%----------------------------
%---begin part 1-------------
%----------------------------

% convert to baseband QAM signal
t = (0:length(xRF)-1)'*Ts; %time of xRF
N = length(cp);

xBB = 2*exp(-1j*2*pi*fc*t).*xRF; % desired baseband signal

% examine spectrum of unfiltered baseband
figure
subplot(1,4,1)
spec_analysis(xBB,fs);
title("Part I - unfiltered baseband signal xBB");

%filter out out-of-spectrum components
pR=pT; % receiver filter
xBBf = conv(xBB, pR); % filtered baseband signal

subplot(1,4,2)
spec_analysis(xBBf,fs);
title("Part I - filtered baseband signal xBBf");
% examine spectrum of filtered baseband

% look at eye pattern
subplot(1,4,3)
eye_pattern(xBBf);
title("Part I - eye pattern of xBBf")

% non-data aided timing recovery (10.1)
% sample the signal in timing phase that maximizes signal power
% find max of rhot for sample locations for decimation
L = Tb/Ts;
rhot = find_rhot(xBBf, L);
[value, indices] = max(rhot);
peak_index = indices(1);
xBBd = xBBf(peak_index:L:end); % decimated baseband signal

% identify preamble and extract payload
preamble = [cp; cp; cp; cp];
payload_start = find_payload_start(xBBd, cp);

payload = xBBd(payload_start:end);
subplot(1,4,4)
eye_pattern(payload)
title("Part I - eye pattern of payload")

% convert QAM siganl to bit array
bits = QPSK2bits(payload); % data bits

% save bits to file
bin2file(bits', "transmitted_file.txt");

%----------------------------
%---end of part 1------------
%----------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Determine Fourier series coeffient rho_n from signal cBB        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pn = find_rhon(cBB, Tb, n)
    sigma_s2 = var(cBB);
    CBB = fft(cBB);
    pn = sigma_s2/Tb*trapz(CBB.*delayseq(CBB, n/Tb)'.'); % Equation 10.7
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Determine time recovery function rho(t) from signal cBB        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rhot = find_rhot(cBB, L)
    rhot = zeros(1);
    for tau = 1:L
        rhot(tau) = sum(abs(cBB(tau:L:end)).^2);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Find start of preamble of signal y with pilot cp                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function index = find_preamble_start(y, cp)
    correlations = zeros(1,length(y));
    for i = 1:1:length(y)
        if i+length(cp)-1 < length(y)
            % find correlation of window of signal y with pilot cp 
            correlations(i) = corr(real(y(i:i+length(cp)-1)), real(cp));
        end
    end

    % find beginning of window with highest correlation
    index = find(correlations>0.9, 1, "first");
end

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

function eye_pattern(y)
    plot(y, "b")
    xlabel("real part")
    ylabel("imaginary part");
end

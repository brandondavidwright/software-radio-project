close all;

load('xRF4.mat')

%-----start of part 1-----

% examine spectrum of xRF
figure
spec_analysis(xRF, fs)
title("Received signal");

% convert to baseband QAM signal
t = (0:length(xRF)-1)'*Ts; %time of xRF
xBB = 2*exp(-1j*2*pi*fc*t).*xRF; % desired baseband signal 

%filter out out-of-spectrum components
pR=pT; % receiver filter
xBBf = conv(xBB, pR);
% construct preamble
preamble = [cp; cp; cp; cp];
N = length(cp);

% Allow for variable timing phase
timing_phase_start = 23;

% Decimate at twice the symbol rate
M = 2;
L = Tb/Ts/M; % decimation factor L
xBBd = decimator(xBBf(timing_phase_start:end),L); % baseband signal decimated at twice the symbol rate

abs_ryy = abs(autocorrelation(xBBd, M*N-1)); % N = 31
figure
plot(abs_ryy)
xlabel("n")
ylabel("autocorrelation")
title("Autocorrelation");

% find peaks for fractionally-space equalizer
ryy_maxima_indeces = find(islocalmax(abs_ryy)==1);
ryy_max_index = find_autocorrelation_peak(abs_ryy, N, M);
s = xBBd(ryy_max_index:ryy_max_index+M*N-1);
% estimate tap weigths
[w, mse] = estimate_tap_weigths(s, cp, N, M);
%center largest tap weight
w = circshift(w, M*N/2 - find(w==max(w)));
% equalize with fractionally-spaced equalizer
xBBe = fractionally_spaced_eq(xBBd, w, M*N); %equalized baseband siganl
figure
eye_pattern(xBBe);
title("xBBe")

payload_start = find_payload_start(xBBe, cp, M);
%decimate payload by 2
payload = decimator(xBBe(payload_start:end),2);
figure
subplot(1,2,1)
eye_pattern(xBBf)
title("Part III - filtered baseband signal of xRF5")
subplot(1,2,2)
eye_pattern(payload);
title("Part III - payload of xRF5")

% convert QAM siganl to bit array
bits = QPSK2bits(payload); % data bits

% save bits to file
bin2file(bits', "transmitted_file.txt");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Find start of preamble of signal y with pilot cp                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function index = find_preamble_start(y, cp, M)
    correlations = zeros(1,length(y));
    for i = 1:1:length(y)
        if i+M*length(cp)-1 < length(y)
            correlations(i) = corr(real(y(i:M:i+M*length(cp)-1)), real(cp));
        end
    end

    index = find(correlations==max(correlations));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Find start of payload of signal y and pilot cp                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function index = find_payload_start(y, cp, M)
    correlations = zeros(1,length(y));
    for i = 1:1:length(y)
        if i+M*length(cp)-1 < length(y)
            % find correlation of wind of signal y with pilot cp
            correlations(i) = corr(real(y(i:M:i+M*length(cp)-1)), real(cp));
        end
    end
    % Find beginning of last correlated pilot segment
    start_last_cp = find(correlations>0.9, 1, "last");
    % Find beginning of payload
    index = start_last_cp+M*length(cp)-1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Determine autocorrelation of signal y                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ryy = autocorrelation(y, N)
   ryy = zeros(1);
   for n = 1:1:length(y)-N+1
      r=0;
      for k = 0:1:N+1
         if n-k > 0 && n-N-1-k>0
            % Find running autocorrelation
            r = y(n-k)*conj(y(n-N-1-k)) + r; % Equation from Problem 10.10c
         end
      end
      ryy(n) = r;      
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Estimate tap weights for equalizer from pilot s in segment y    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [w, e] = estimate_tap_weigths(y, s, N, M)
    yi = fliplr(y);
    mu = 0.001; %lower mu due to non symbols in tap weight estimation
    iterations = M*100000;
    w = zeros(M*N,1);
    e = zeros(iterations+1,1);
    % NLMS adaptation algorithm
    for i = 0:iterations
        ei = s(floor(mod(i, M*N)/M)+1) - w'*yi;
        e(i+1) = ei;
        w = w + 2*mu*conj(ei)*yi;
        yi = circshift(yi,-1);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Fractionally-spaced equalizer of signal y of order N            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s2 = fractionally_spaced_eq(y, w, N)
    overshoot = N-(mod(length(y),N)+1);
    %equalize signal based on tap weights
    for n = 1:1:length(y)
        if n + N-1 > length(y)
            s2(n:length(y)) = w(1:mod(length(y),N)+overshoot)'*y(n:length(y));
            overshoot = overshoot - 1;
        else
            s2(n:n+N-1) = w'*y(n:n+N-1);
        end
    end
    s2 = reshape(s2,length(s2),1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Determine peak in autocorrelation of signal y                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function peak_index = find_autocorrelation_peak(abs_ryy, N, M)
    L = N*M;
    before_payload = find(abs_ryy>0);
    begin_preamble = before_payload(1)+1;
    peak = max(abs_ryy(begin_preamble:begin_preamble+2*L+N));
    peak_index = find(abs_ryy(1:begin_preamble+2*L+N)==peak);
end

function eye_pattern(y)
    plot(y, "b")
    xlabel("real part")
    ylabel("imaginary part");
end

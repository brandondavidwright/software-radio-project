close all; clear;

load('xRF5.mat')

%----------------------------
%---begin part 1-------------
%----------------------------

% examine spectrum of xRF
% figure
% spec_analysis(xRF, fs)
% title("Received signal");

% convert to baseband QAM signal
t = (0:length(xRF)-1)'*Ts; %time of xRF
N = length(cp);

xBB = 2*exp(-1j*2*pi*fc*t).*xRF; % desired baseband signal

% examine spectrum of unfiltered baseband
% figure
% spec_analysis(xBB,fs);
% title("Unfiltered baseband signal");

%filter out out-of-spectrum components
pR=pT; % receiver filter
xBBf = conv(xBB, pR); % filtered baseband signal

% figure
% spec_analysis(xBBf,fs);
% title("filtered baseband signal");
% examine spectrum of filtered baseband

% non-data aided timing recovery (10.1)
% sample the signal in timing phase that maximizes signal power
% sample the signal in timing phase that maximizes signal power
% find Fourier series coeffients for rho(n)
rho0 = find_rhon(xBBf, Tb, 0);
rho1 = find_rhon(xBBf, Tb, 1);

% find timing cost recovery function rho(t)
t1 = (0:length(xBBf)-1)'*Ts;
rhot = rho0 + 2*abs(rho1)*cos(2*pi/Tb*t1 + angle(rho1)); % 10.12

%find max power of periodic rho(t)
peak_indices = find(rhot==max(rhot));
L = Tb/Ts;
xBBd = xBBf(peak_indices); % decimated baseband signal %

% figure
% eye_pattern(xBBd)
% title("eye pattern of xBBd")

% identify preamble and extract payload
preamble = [cp; cp; cp; cp];
N = length(cp);

%----------------------------
%---end of part 1------------
%----------------------------
%---start of part 2----------
% equalize channel
% determine autocorrelation to find beginning of pilot sequence
abs_ryy = abs(autocorrelation(xBBd, N-1)); % N = 31
figure
plot(abs_ryy)
xlabel("n")
ylabel("|r_y_y|")
title("Part II - autocorrelation of xRF2");

% locate peak to find pilot sequence
peak_index = find_starting_peak(abs_ryy, N);
s = xBBd(peak_index:peak_index+N-1); % pilot
% estimate tap weights for equalizer
w = estimate_tap_weigths(s, cp, N);
% put highest tap weight in center
w = circshift(w, N/2 - find(w==max(w)));
% equalize with symbol-spaced equalizer
xBBe = symbol_spaced_equalizer(xBBd, w, N); % equalized baseband signal

% figure
% eye_pattern(xBBe);
% title("xBBe");

% find begining of payload
payload_start = find_payload_start(xBBe, cp);
% extract payload
payload_p2 = xBBe(payload_start:end);

% look at eye pattern
figure
subplot(1,2,1)
eye_pattern(xBBf);
title("Part II - eye pattern of filtered baseband for xRF5")
subplot(1,2,2)
eye_pattern(payload_p2);
title("Part II - eye pattern of payload of xRF5")

% convert QAM siganl to bit array
bits = QPSK2bits(payload_p2); % data bits

% save bits to file
bin2file(bits', "transmitted_file.txt");

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
        rhot(tau) = sum(abs(cBB(tau:1:end)).^2);
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
function w = estimate_tap_weigths(y, s, N)
    yi = fliplr(y);
    mu = 0.0025;
    iterations = 100000;
    w = zeros(32,1);
    
    % NLMS adaptation algorithm ***TODO is this right?
    for i = 1:1:iterations
        ei = s(mod(i, N)+1) - w'*yi;
        w = w + 2*mu*conj(ei)*yi;
        yi = circshift(yi,-1);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Sybmol-spaced equalizer of signal y with tap weights w           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s2 = symbol_spaced_equalizer(y, w, N)
    % find remainder of length of y % N
    remainder = N-(mod(length(y),N)+1);
    %equalize signal based on tap weights
    for n = 1:1:length(y)
        if n + N-1 > length(y)
            s2(n:length(y)) = w(1:mod(length(y),N)+remainder)'*y(n:length(y));
            remainder = remainder - 1;
        else
            s2(n:n+N-1) = w'*y(n:n+N-1);
        end
    end
    s2 = reshape(s2,length(s2),1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Find peak in autocorrelation to find segment corresponding to cp%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function index = find_starting_peak(y, N)
    before_payload = find(y>0);
    begin_preamble = before_payload(1)+1;
    index = find(y==max(y(begin_preamble:begin_preamble+2*N)));
end

function eye_pattern(y)
    plot(y, "b")
    xlabel("real part")
    ylabel("imaginary part");
end

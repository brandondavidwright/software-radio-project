close all; clear;

load('xRF5.mat')

%-----start of part 1-----

% examine spectrum of xRF
figure
spec_analysis(xRF, fs)
title("Received signal");

% convert to baseband QAM signal
t = (0:length(xRF)-1)'*Ts; %time of xRF
xBB = 2*exp(-1j*2*pi*fc*t).*xRF; % desired baseband signal 

% examine spectrum of unfiltered baseband
% figure
% spec_analysis(xBB,fs);
% title("Unfiltered baseband signal");

%filter out out-of-spectrum components
pR=pT; % receiver filter
xBBf = conv(xBB, pR);
% construct preamble
preamble = [cp; cp; cp; cp];
N = length(cp);

% examine spectrum of unfiltered baseband
% figure
% spec_analysis(xBBf,fs);
% title("Filtered baseband signal");

% Decimate at twice the symbol rate
M = 2;
L = Tb/Ts/M; % decimation factor L
xBBd = decimator(xBBf,L); % baseband signal decimated at twice the symbol rate
% ***how do we vary timing phase? TODO***

% look at eye pattern
% figure
% eye_pattern(xBBd);
% title("Eye pattern of xBBd")
%---start of part 2----------
% equalize channel
% determine autocorrelation to find beginning of pilot sequence

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
w = estimate_tap_weigths(s, cp, N, M);
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
function w = estimate_tap_weigths(y, s, N, M)
    yi = fliplr(y);
    mu = 0.001; %lower mu due to non symbols in tap weight estimation
    iterations = M*100000;
    w = zeros(M*N,1);
    
    % NLMS adaptation algorithm ***TODO is this right?
    for i = 0:iterations
        ei = s(floor(mod(i, M*N)/M)+1) - w'*yi;
        w = w + 2*mu*conj(ei)*yi;
        yi = circshift(yi,-1);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Fractionally-spaced equalizer of signal y of order N            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s2 = fractionally_spaced_eq(y, w, N)
    overshoot = N-(mod(length(y),N)+1); %calculate this
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
%  Determine peak in autocorrelation of signal y                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function peak_index = find_autocorrelation_peak(abs_ryy, N, M)
    L = N*M;
    before_payload = find(abs_ryy>0);
    begin_preamble = before_payload(1)+1;
    peak = max(abs_ryy(begin_preamble:begin_preamble+2*L+N));
    peak_index = find(abs_ryy(1:begin_preamble+2*L+N)==peak);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Estimate phase offset phic of signal y                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function phic = esimate_phase_offset(y)
    phi = zeros(size(y));
    s1 = zeros(size(y));
    mu = 0.001;

    % decision-directed phase recovery loop
    for n = 1:1:length(y)-1
        s1(n) = y(n)*exp(-1j*phi(n));
        s2 = sign(real(s1(n)))+1j*sign(imag(s1(n))); % slicer
        s12=s1(n)*s2';
        e = imag(s12)/real(s12);  % Equation 9.50
        phi(n+1) = phi(n) + mu*e;
    end
    phic = phi(end);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Determine complex sum J of signal y for interval N1 to N2       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function J = findJ(y, N1, N2, N, M)
    J = 0;
    for n = N1:1:N2
        J = J + y(n + M*N)*y(n)';
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Find end of transicence in signal y with period N               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- https://www.mathworks.com/matlabcentral/fileexchange/68486-transient-time-measurement-with-matlab-implementation
% --- TODO: cite this
function index = find_end_transience(y, N)
    for i = 1:1:150
        dev = abs(abs(y(i)) - abs(y(i+N)))/abs(y(i));
        if dev < 0.05 %got 0.02 value from link above, but 0.05 seems to work better for both test files
            break;
        end
    end
    index = i;
end

function eye_pattern(y)
    plot(y, "b")
    xlabel("real part")
    ylabel("imaginary part");
end

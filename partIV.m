close all; clear;

load('xRF8.mat')

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
% find Fourier series coeffients for rho(n)
rho0 = find_rhon(xBBf, Tb, 0);
rho1 = find_rhon(xBBf, Tb, 1);

% find timing cost recovery function rho(t)
t1 = (0:length(xBBf)-1)'*Ts;
rhot = rho0 + 2*abs(rho1)*cos(2*pi/Tb*t1 + angle(rho1)); % 10.12

%find max power of periodic rho(t)
peak_indices = find(rhot==max(rhot));

xBBd = xBBf(peak_indices); % decimated baseband signal 

figure
eye_pattern(xBBd)
title("eye pattern of xBBd")
%----------------------------
%---part 1 to be continued---
%----------------------------
%---begin part 4.1-----------
%----------------------------
% estimate delta_fc with
% pilot-aided carrier recovery (9.3)
% find interval N1 to N2 after initial transient interval
N1 = find_end_transience(xBBd, N);
N2 = N1 + N - 1;
% find J
J = findJ(xBBd, N1, N2, N);
% estimate delta_fc
Dfc_est = angle(J)/(2*pi*N*Tb); % equation 9.47
t1 = (0:1:length(xBBd)-1)'*Tb;
% adjust carrier frequency with estimate
offset = exp(-1j*2*pi*Dfc_est*t1);
xBBo = xBBd.*offset; % baseband signal adjusted for carrier offset

figure
eye_pattern(xBBo);
title("eye pattern of xBBo");

%----------------------------
%---end part 4.1-------------
%----------------------------
%---part 1 continues---------
%----------------------------

% identify preamble and extract payload
preamble = [cp; cp; cp; cp];
N = length(cp);

%----------------------------
%---end of part 1------------
%----------------------------
%---start of part 2----------
% equalize channel
% determine autocorrelation to find beginning of pilot sequence
abs_ryy = abs(autocorrelation(xBBo, N-1)); % N = 31
figure
plot(abs_ryy)
xlabel("n")
ylabel("autocorrelation")
title("Autocorrelation");

% locate peak to find pilot sequence
peak_index = find_starting_peak(abs_ryy, N);
s = xBBo(peak_index:peak_index+N-1); % pilot
% estimate tap weights for equalizer
w = estimate_tap_weigths(s, cp, N);
% put highest tap weight in center
w = circshift(w, N/2 - find(w==max(w)));
% equalize with symbol-spaced equalizer
xBBe = symbol_spaced_equalizer(xBBo, w, N); % equalized baseband signal

figure
eye_pattern(xBBe);
title("xBBe");
%----------------------------
%---part 2 to be continued---
%----------------------------
%---begin part 4.3-----------
%----------------------------
% estimate phase offset phi_c
phic = esimate_phase_offset(xBBe);

xBBp = xBBe*exp(-1j*phic); % phase-adjusted baseband signal

figure
eye_pattern(xBBp);
title("Eye view of xBBp")

%----------------------------
%---end part 4.3-------------
%----------------------------
%---part 2 continues---------

% find begining of payload
payload_start = find_payload_start(xBBp, cp);
% extract payload
payload = xBBp(payload_start:end);

% look at eye pattern
figure
subplot(1,2,1)
eye_pattern(xBBf);
title("Part IV - filtered baseband of xRF8")
subplot(1,2,2)
eye_pattern(payload);
title("Part IV - payload of xRF8")

% convert QAM siganl to bit array
bits = QPSK2bits(payload); % data bits

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
    mu = 0.001;
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
function J = findJ(y, N1, N2, N)
    J = 0;
    for n = N1:1:N2
        J = J + y(n + N)*y(n)';
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Find end of transicence in signal y with period N               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- https://www.mathworks.com/matlabcentral/fileexchange/68486-transient-time-measurement-with-matlab-implementation
% --- TODO: cite this
function index = find_end_transience(y, N)
    for i = 1:1:150
        dev = abs(abs(y(i)) - abs(y(i+N)));
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

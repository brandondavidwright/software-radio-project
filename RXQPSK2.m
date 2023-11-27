close all; clear;

load('xRF1.mat')

%-----start of part 1-----

% examine spectrum of xRF
figure
spec_analysis(xRF, fs)
title("Received signal");

% convert to baseband QAM signal
t = (0:length(xRF)-1)'*Ts; %time of xRF
xBB = 2*exp(-1j*2*pi*fc*t).*xRF; % desired baseband signal

% examine spectrum of unfiltered baseband
figure
spec_analysis(xBB,fs);
title("Unfiltered baseband signal");

%filter out out-of-spectrum components
pR=pT; % receiver filter
y = conv(xBB, pR);

% examine spectrum of unfiltered baseband
figure
spec_analysis(y,fs);
title("Filtered baseband signal");

%--remove timing phase for part 3
% % sample the signal in timing phase that maximizes signal power
% % 10.1
% [CBB, f] = spec_analysis(y,fs);
% 
% rho0 = 0;
% % ***is this the right way to delay?
% rho1 = find_rho1(CBB, Tb);
% 
% % calculate timing recovery cost function
% rhot = rho0 + 2*abs(rho1)*cos(2*pi/Tb*t + angle(rho1)); % 10.12
% sample_indeces = find(rhot==max(rhot));
% 
% % figure this part out
% xBBd = y(sample_indeces); % max power

%------start part 3------
M = 2;
L = 100/M;
xBBd = decimator(y,L); % decimated at twice the symbol rate
% ***how do we vary timing phase?

% look at eye pattern
figure
eye_pattern(y);

% identify preamble and extract payload

%detect preamble
% extract data symbols from payload and convert to bit stream
preamble = [cp; cp; cp; cp];
N = length(cp);
%payload_start = find_payload_start(xBBd, cp);

%payload_p1 = xBBd(payload_start:end);
%figure
%eye_pattern(payload_p1)

%bits = QPSK2bits(payload);

%-------end of part 1--------
%----------------------------
%-------start of part 2------
abs_ryy = abs(autocorrelation(xBBd, M*(N-1))); % N = 31
abs_ryy_dec = expander(abs(autocorrelation(decimator(xBBd(1:end),M), N-1)),M);
figure
plot(abs_ryy)
xlabel("n")
ylabel("autocorrelation")
title("Autocorrelation");

ryy_maxima_indeces = find(islocalmax(abs_ryy)==1);
ryy_max_index = find_autocorrelation_peak(abs_ryy, N, M); %should be 199 for RF2, 149 for RF1
s = xBBd(ryy_max_index:ryy_max_index+M*N-1);
w = esimte_tap_weights(s, cp, N, M);
%center largest tap weight
w = circshift(w, N/2 - find(w==max(w)));

xBBe = fractionally_spaced_eq(xBBd, w, M*N);
figure
eye_pattern(xBBe);
title("xBBe")

N = length(cp);
payload_start = find_payload_start(xBBe, cp, M);

payload_p2 = decimator(xBBe(payload_start:end),2);
figure
eye_pattern(payload_p2);
title("payload part of xBBe")

bits = QPSK2bits(payload_p2); % TODO fix this

% save bits to file
bin2file(bits', "transmitted_file.txt");

function p1 = find_rho1(CBB, Tb)
    sigma_s = var(CBB); 
    p1 = 0;
    for f = 1:1:length(CBB)
        if(f>1/Tb+1)
            p1 = p1 + sigma_s^2/Tb*CBB(f)*CBB(f-1/Tb)';
        end
    end
end

function index = find_preamble_start(y, cp, M)
    for i = 1:1:length(y)
        if i+M*length(cp)-1 < length(y)
            correlations(i) = corr(real(y(i:M:i+M*length(cp)-1)), real(cp));
        end
    end

    index = find(correlations==max(correlations));
end

function index = find_payload_start(y, cp, M)
    for i = 1:1:length(y)
        if i+M*length(cp)-1 < length(y)
            correlations(i) = corr(real(y(i:M:i+M*length(cp)-1)), real(cp));
        end
    end
    start_last_cp = find(correlations>0.9, 1, "last");
    
    index = start_last_cp+M*length(cp)-1;
end

function ryy = autocorrelation(y, N)
   ryy = zeros(1);
   for n = 1:1:length(y)-N+1
      r=0;
      for k = 0:1:N+1
         if n-k > 0 && n-N-1-k>0
            r = y(n-k)*conj(y(n-N-1-k)) + r;
         end
      end
      ryy(n) = r;      
    end
end

function w = esimte_tap_weights(y, s, N, M)
    % yi = flipud(y); %TA said to flip was he wrong?
    yi = y;
    mu = 0.001; %lower mu due to non symbols in tap weight estimation
    iterations = M*100000;
    w = zeros(M*N,1);
    
    for i = 0:iterations        
        ei = s(floor(mod(i, M*N)/M)+1) - w'*yi;
        w = w + 2*mu*conj(ei)*yi;
        yi = circshift(yi,-1);
    end
end

function s2 = fractionally_spaced_eq(y, w, N)
    % add adaptation algorithm?  Like equalizerT2_NLMS.m',
    w_row = w;
    overshoot = N-(mod(length(y),N)+1); %calculate this
    for n = 1:1:length(y)
        if n + N-1 > length(y)
            s2(n:length(y)) = w_row(1:mod(length(y),N)+overshoot)'*y(n:length(y));
            overshoot = overshoot - 1;
        else
            s2(n:n+N-1) = w_row'*y(n:n+N-1);
        end
    end
    s2 = reshape(s2,length(s2),1);
end

function s2 = symbol_spaced_eq(y, w, N)
    w_row = w;
    overshoot = N -(mod(length(y),N)+1);
    for n = 1:1:length(y)
        if n + N-1 > length(y)
            s2(n:length(y)) = w_row(1:mod(length(y),N)+overshoot)'*y(n:length(y));
            overshoot = overshoot - 1;
        else
            s2(n:n+N-1) = w_row'*y(n:n+N-1);
        end
    end
    s2 = reshape(s2,length(s2),1);
end

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

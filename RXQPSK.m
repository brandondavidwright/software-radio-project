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

% sample the signal in timing phase that maximizes signal power
% 10.1
[CBB, f] = spec_analysis(y,fs);

rho0 = 0;
% ***is this the right way to delay?
rho1 = find_rho1(CBB, Tb);

% calculate timing recovery cost function
rhot = rho0 + 2*abs(rho1)*cos(2*pi/Tb*t + angle(rho1)); % 10.12
sample_indeces = find(rhot==max(rhot));

% figure this part out
xBBd = y(sample_indeces); % max power

% look at eye pattern
figure
eye_pattern(y);

% identify preamble and extract payload

%detect preamble
% extract data symbols from payload and convert to bit stream
preamble = [cp; cp; cp; cp];
N = length(cp);
%preamble_end = find_preamble_start(xBBd, preamble);
payload_start = find_payload_start(xBBd, cp);

%payload = xBBd(preamble_end + length(preamble):end);
payload = xBBd(payload_start:end);
figure
eye_pattern(payload)

%bits = QPSK2bits(payload);

%-------end of part 1--------
%----------------------------
%-------start of part 2------
abs_ryy = abs(autocorrelation(xBBd, N-1)); % N = 31
figure
plot(abs_ryy)
xlabel("n")
ylabel("autocorrelation")
title("Autocorrelation");

ryy_maxima_indeces = find(islocalmax(abs_ryy)==1);
%**** TODO - which peak should I pick? *****
% does cp = s?
s = xBBd(ryy_maxima_indeces(4):ryy_maxima_indeces(4)+N-1);
w = ssce(s, cp, N);
%center largest tap weight
w = circshift(w, N/2 - find(w==max(w)));

xBBe = equalizer(xBBd, w, N);
figure
eye_pattern(xBBe);
title("xBBe")

N = length(cp);
payload_start = find_payload_start(xBBe, cp);

payload = xBBe(payload_start:end);
figure
eye_pattern(payload);
title("payload part of xBBe")

bits = QPSK2bits(payload); % TODO fix this

% save bits to file
bin2file(bits', "transmitted_file.txt");

function p1 =  find_rho1(CBB, Tb)
    sigma_s = var(CBB); 
    p1 = 0;
    for f = 1:1:length(CBB)
        if(f>1/Tb+1)
            p1 = p1 + sigma_s^2/Tb*CBB(f)*CBB(f-1/Tb)';
        end
    end
end

function index = find_preamble_start(y, preamble)
    for i = 1:1:length(y)
        if i+length(preamble)-1 < length(y)
            correlations(i) = corr(real(y(i:i+length(preamble)-1)), real(preamble));
        end
    end

    index = find(correlations==max(correlations));
end

function index = find_payload_start(y, cp)
    for i = 1:1:length(y)
        if i+length(cp)-1 < length(y)
            correlations(i) = corr(real(y(i:i+length(cp)-1)), real(cp));
        end
    end
    start_last_cp = find(correlations>0.9, 1, "last");
    
    index = start_last_cp+length(cp);
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

function w = ssce(y, s, N)
    yi = fliplr(y);
    mu = 0.005;
    iterations = 1000000;
    w = zeros(32,1);
    
    for i = 1:1:iterations
        ei = s(mod(i, N)+1) - w'*yi;
        w = w + 2*mu*conj(ei)*yi;
        yi = circshift(yi,-1);
    end
end

function s2 = equalizer(y, w, N)
    s2 = zeros(1);
    for n = 1:1:length(y)
        if n + N-1 > length(y)
            s2(n:length(y)) = y(n:length(y));
        else
            s2(n:n+N-1) = w'*y(n:n+N-1);
        end
    end
    s2 = s2';
end

function eye_pattern(y)
    plot(y, "b")
    xlabel("real part")
    ylabel("imaginary part");
end

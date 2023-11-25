close all; clear;

load('xRF1.mat')

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
plot(y, "b")
xlabel("real part")
ylabel("imaginary part");

% identify preamble and extract payload

%detect preamble
% extract data symbols from payload and convert to bit stream
preamble = [cp; cp; cp; cp];

preamble_start = find_preamble_start(xBBd, preamble);

payload = xBBd(preamble_start + length(preamble):end);

bits = QPSK2bits(payload);

% save bits to file
bin2file(bits, "transmitted_file.txt");

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
    correlations = zeros(1);

    for i = 1:1:length(y)
        if i+length(preamble)-1 < length(y)
            correlations(i) = corr(y(i:i+length(preamble)-1), preamble);
        end
    end

    index = find(correlations==max(correlations));
end

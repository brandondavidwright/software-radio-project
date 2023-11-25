close all; clear;

load('xRF1.mat')

% convert to baseband QAM signal
xT = xRF;

% QAM demodulation
t = (0:length(xRF)-1)'*Ts; %time of xRF
xBB = 2*exp(-1j*2*pi*fc*t).*xT; % desired baseband signal

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
CBB = fft(y);

sigma_s = 1;
rho0 = 0;
% ***is this the right way to delay?
rho1 = sigma_s/Tb*cumtrapz(CBB*delayseq(1/Tb, CBB')); % 10.9

% calculate timing recovery cost function
rhot = rho0 + 2*abs(rho1)*cos(2*pi/Tb*t + angle(rho1)); % 10.12

% figure this part out
xBBd = 0; % max power

% look at eye pattern
figure
plot(y, "b")
xlabel("real part")
ylabel("imaginary part");

% identify preamble and extract payload

% extract data symbols from payload and convert to bit stream
N=32;
payload = xBBd(4*32:end);
b = QPSK2bits(xBBd);

% save bits to file
bin2file(b, "transmitted_file.gif");

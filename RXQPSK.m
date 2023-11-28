close all; clear;

load('xRF8.mat')

%-----start of part 1-----

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
xR = conv(xBB, pR);

figure
spec_analysis(xR,fs);
title("filtered baseband signal");
% examine spectrum of unfiltered baseband

% look at eye pattern
figure
eye_pattern(xR);
title("eye pattern of xR")

% sample the signal in timing phase that maximizes signal power
% 10.1
CBB = fft(xR);

rho0 = find_rhon(CBB, Tb, 0);
rho1 = find_rhon(CBB, Tb, 1);
% 
% % calculate timing recovery cost function
rhot = rho0 + 2*abs(rho1)*cos(2*pi/Tb*t + angle(rho1)); % 10.12
sample_indices = find(rhot==max(rhot));

% figure this part out
xBBd = xR(sample_indices); % max power

figure
eye_pattern(xBBd)
title("eye pattern of xBBd")
%------------------
%---begin part 4.1-
%------------------
i = find_end_transience(xBBd, N);
N1=find_end_transience(xBBd, N)-2; %18 for xrf7 - why subtract 2?
N2 = N1 + N -1;
J = findJ(xBBd, N1, N2, N);

Dfc_est = angle(J)/(2*pi*N*Tb);
t1 = (0:1:length(xBBd)-1)'*Tb;
offset = exp(-1j*2*pi*Dfc_est*t1);
xBBo = xBBd.*offset;
figure
eye_pattern(xBBo);
title("eye pattern of xBBo");

%------------------
%---end part 4.1---
%------------------

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
abs_ryy = abs(autocorrelation(xBBo, N-1)); % N = 31
figure
plot(abs_ryy)
xlabel("n")
ylabel("autocorrelation")
title("Autocorrelation");

peak_index = find_starting_peak(abs_ryy, N);
s = xBBo(peak_index:peak_index+N-1);
w = estimate_tap_weigths(s, cp, N);
w = circshift(w, N/2 - find(w==max(w)));

xBBe = symbol_spaced_equalizer(xBBo, w, N);
figure
eye_pattern(xBBe);
title("xBBe");


%------------------
%--begin part 4.3---
%------------------

phic = ddrc(xBBe);

xBBc = xBBe*exp(-1j*phic); %TODO negative or positive?

figure
eye_pattern(xBBc);
title("Eye view of xBBc")

%------------------
%--end part 4.3-----
%------------------

payload_start = find_payload_start(xBBc, cp);

payload_p2 = xBBc(payload_start:end);
figure
eye_pattern(payload_p2);
title("payload")

bits = QPSK2bits(payload_p2); % TODO fix this

% save bits to file
bin2file(bits', "transmitted_file.txt");

function pn = find_rhon(CBB, Tb, n)
    sigma_s = var(CBB); 
    pn = 0;
    % for f = 1:1:length(CBB)
    %     if(f>1/Tb+1)
    %         pn = pn + sigma_s^2/Tb*CBB(f)*CBB(f-n/Tb)';
    %     end
    % end
    pn = trapz(CBB.*reshape(delayseq(CBB, n/Tb)',length(CBB),1));
end

function index = find_preamble_start(y, cp)
    correlations = zeros(1,length(y));
    for i = 1:1:length(y)
        if i+length(cp)-1 < length(y)
            correlations(i) = corr(real(y(i:i+length(cp)-1)), real(cp));
        end
    end

    index = find(correlations>0.9, 1, "first");
end

function index = find_payload_start(y, cp)
    correlations = zeros(1,length(y));
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

function w = estimate_tap_weigths(y, s, N)
    % yi = flipud(y); %TA said to flip
    yi = y;
    mu = 0.0025;
    iterations = 100000;
    w = zeros(32,1);
    
    for i = 1:1:iterations
        ei = s(mod(i, N)+1) - w'*yi;
        w = w + 2*mu*conj(ei)*yi;
        yi = circshift(yi,-1);
    end
end

function s2 = symbol_spaced_equalizer(y, w, N)
    overshoot = N-(mod(length(y),N)+1); %calculate this
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

function dfc = ddrc(y)
    phi = zeros(size(y));
    s1 = zeros(size(y));
    mu = 0.001;
    for n = 1:1:length(y)-1
        s1(n) = y(n)*exp(-1j*phi(n));
        s2 = sign(real(s1(n)))+1j*sign(imag(s1(n))); %slicer
        s12=s1(n)*s2';
        e = imag(s12)/real(s12);
        phi(n+1) = phi(n) + mu*e;
    end
    dfc = phi(end);
end

function J = findJ(y, N1, N2, N)
    J = 0;
    for n = N1:1:N2
        J = J + y(n + N)*y(n)';
    end
end

% --- https://www.mathworks.com/matlabcentral/fileexchange/68486-transient-time-measurement-with-matlab-implementation
% --- TODO: cite this
function index = find_end_transience(y, N)
    for i = 1:1:150
        dev = abs(abs(y(i)) - abs(y(i+N)));
        disp(dev)
        disp(i)
        if dev < 0.02 %got this value from link above
            index = i;
            break;
        end
    end
end

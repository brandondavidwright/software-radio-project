alpha=0.25;N=36;M=6;
p_T=sr_cos_p(N,M,alpha)';                   %Transmit filter
p_R=p_T;                                    %Receive filter
fs=(1/M)*(1+alpha);
c=firpm(100,[0 fs 1.25*fs 1],[1 1 0 0])';   %Channel
h=conv(c,conv(p_T,p_R));                    %0verall response
h=h/max(h);                                 %Normalized to the peak of 1 
[x,i]=max(h);                               %Take the relevant samples of
h=h(i-floor(length(h)/M/2)*M:M:end);        %The overall response
b=sign(randn(60000,1));
for k=1:length(b)/6;                        %Generate 64-QAM symbols
    i=(k-1)*6;
    s(k)=4*b(i+1)+2*b(i+2)+b(i+3)+j*(4*b(i+4)+2*b(i+5)+b(i+6));
end
x=filter(h,1,s);                            %Pass data symbols through the cahnnel
x=x(round(length(h)):end);
subplot(1,2,1),plot(x,'.')
axis('square')
axis([-8 8 -8 8])
xlabel('(a)')
subplot(1,2,2),plot(s,'.')
axis('square')
axis([-8 8 -8 8])
xlabel('(b)')
%axes('position',[0.25 0.25 0.5 0.5]);




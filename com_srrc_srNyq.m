
N=29;M=5;alpha=0.25;gamma=0.1;
hsr=sr_cos_p(N,M,alpha);
f2=(1/M)*(1+alpha);
hlp=firpm(100,[0 f2 1.25*f2 1],[1 1 0 0])';
h1=conv(hsr,hsr);
h1=conv(h1,hlp);
hsrNyq=sr_Nyquist_p(N,M,alpha,gamma);
h2=conv(hsrNyq,hsrNyq);
h2=conv(h2,hlp);
h1=h1/max(h1);
h2=h2/max(h2);
[x,i]=max(h1);
h1=h1(i-floor(length(h1)/(2*M))*M:M:end)
[x,i]=max(h2);
h2=h2(i-floor(length(h2)/(2*M))*M:M:end)

b=sign(randn(60000,1));
for k=1:length(b)/6;
    i=(k-1)*6;
    s(k)=4*b(i+1)+2*b(i+2)+b(i+3)+j*(4*b(i+4)+2*b(i+5)+b(i+6));
end
x1=filter(h1,1,s);
x2=filter(h2,1,s);

Hsr=20*log10(abs(fft(hsr,2048)));
Hsr=Hsr(1:1024);
f=[0:1023]/2048;
HsrNyq=20*log10(abs(fft(hsrNyq,2048)));
HsrNyq=HsrNyq(1:1024);
figure(1),axes('position',[0.25 0.25 0.5 0.5])
plot(f,Hsr,'--k',f,HsrNyq,'-k')
legend('square-root raised cosine','square-root Nyquist (M)')
axis([0 0.5 -40 10])
print -deps F4_magresponses.eps
figure(2),axes('position',[0.25 0.25 0.5 0.5])
plot(x1(length(h1):end),'.'),axis('square',[-8 8 -8 8])
print -deps F4_constellation1.eps
figure(3),axes('position',[0.25 0.25 0.5 0.5])
plot(x2(length(h1):end),'.'),axis('square',[-8 8 -8 8])
print -deps F4_constellation2.eps

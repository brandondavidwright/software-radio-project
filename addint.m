function y=addint(x,TsRF,fb,n)
Lx=length(x);
t=[0:Lx-1]*TsRF;
sigmax=std(x);
y=x+5*x.*cos(2*pi*2*fb*t)+4*x.*cos(2*pi*3*fb*t)+n*std(x)*randn(size(x));

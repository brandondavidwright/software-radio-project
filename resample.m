%
% This function resample a signal at a rate slightly different from the
% input rate using linear interpolation between successive samples.
%
% x is the input signal
% epsilon = (fs_out-f_sin)/fs_in, where fs_in is the sampling rate at the
% input and fs_out is the sampling rate at the output.
%
function y=resample(x,epsilon)
t=1;
Dt=1/(1+epsilon);
n=1;n1=1;
y=zeros(round(size(x)*(1+epsilon)));
while n1<length(x)-1
    t=t+Dt;
    n=n+1;
    n1=floor(t);
    n2=ceil(t);
    y(n)=x(n1)+(x(n2)-x(n1))*(t-n1);
end
y=y(1:n);
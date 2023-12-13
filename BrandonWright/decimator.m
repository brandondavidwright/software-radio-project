%
% DECIMATOR: y=decimatorr(x,L) 
% When x is a vector, this function picks one out every L element of x. 
% When x is a matrix, each column of it treated as vector and decimated L
% fold.
function y=decimator(x,L)
[M,N]=size(x);
if (N==1)|(M==1)
    y=x(1:L:end);
else
    y=x(1:L:end,:);
end

% Phase detector with phase unwrapping.
% Variables: 
%   a and b the input variable whose phase difference has to be found.
%   e1 is the previous phase error 
%   e is output phase error unwrapped to be the closest phase anlgle to e1
%
function e=phasedetector(a,b,e1)
e=angle(a*b'); % wrapped phase error 
e=e+2*pi*round((e1-e)/(2*pi)); % add multiples of 2*pi to unwrap the phase
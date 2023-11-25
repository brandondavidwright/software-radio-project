close all; clear;

%
% source
%
info_bits=file2bin('testfile.gif')';      % convert source to bits
%
% generate transmit signal
%     
s=bits2QPSK(info_bits);

b=QPSK2bits(s);

bin2file(b, "bit_test_file.gif");

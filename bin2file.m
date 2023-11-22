% function [data] = bin2file(binData , fname)
% Input(s)
%    fname: is the name of the generated file like 'rec_otello.txt' or '2.bmp'
%    binData: is a column vector bit-stream of data in binary 
%             Sequence of the bytes should be preserved
% Output(s)
%    data: is the what is also written in the file
% Description
%    The function builds the file with the name 'fname' given binary
%    data stream in binData. binData should start with 32 bits representing
%    the length of the file
% Considerations
%    data is zero padded at the end if binData is not of a length of
%    multiple of 8 

function [data] = bin2file(binData , fname)

fid = fopen(fname,'w');
if (mod(length(binData),8) ~= 0) %correct possible missalignment
    binData = [binData ; zeros(8 - mod(length(binData),8),1)]; %padding zeros to the tail
end
mat = reshape(binData,8,length(binData)/8).'; %matrix of the file in binary
mat = num2str(mat);
data = bin2dec(mat); %The byte stream
%extracting the file length from 4 bytes at the beginning
file_len = double(typecast(uint8(data(1:4)),'uint32'));
%write the trimmed data to the file
data = data(5:min(file_len+4,length(data)));
fwrite(fid,data); 
fclose(fid);
return
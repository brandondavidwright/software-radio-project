% function [q] = quan(a , M)
% Input(s)
%    a: the complex input to be quantized
%    M: constellation size
%             
% Output(s)
%    q: decisided symbol
% Description
%   Given the rectangular constellation size , and a point the nearest is
%   in the constellation is returned
% Considerations
%   M must be power of 4 and greater than or equal to 4
%   The constellation coordinates are odd-numbered

function [q] = quan(a , M)

if M == 4
        real_a = real(a);
        imag_a = imag(a);
        if real_a > 0
            if imag_a > 0
                q = 1 + i;
            else
                q = 1 - i;
            end
        else
            if imag_a > 0
                q = -1 + i;
            else
                q = -1 - i;
            end
        end
else
        real_a = real(a);
        imag_a = imag(a);
        if real_a > 0
            if imag_a > 0
                center = (1 + i) * log2(M / 4);
            else
                center = (1 - i) * log2(M / 4);
            end
        else
            if imag_a > 0
                center = (-1 + i) * log2(M / 4);
            else
                center = (-1 - i) * log2(M / 4);
            end
        end
        q = center + quan(a - center, M / 4);
end

return

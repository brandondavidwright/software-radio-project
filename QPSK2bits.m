function info_bits = QPSK2bits(xqpsk)
    L = length(xqpsk);

    info_bits = zeros(1,L);
    
    for k = 1:1:L
        bits = grey_map_decode(xqpsk(k));
        info_bits(2*k - 1) = bits(1);
        info_bits(2*k) = bits(2);
    end
end

function y = grey_map_decode(x) % note: this is different in new version of book
    b00 = -1 + -1j;
    b01 = -1 + 1j;
    b11 = 1 + 1j;
    b10 = 1 - 1j;
    x_round = round(x);
    y = zeros(2,1);
    if(x_round == b00)
        y = [0;0];
    elseif(x_round == b01)
        y = [0;1];
    elseif(x_round == b11)
        y = [1;1];
    elseif(x_round == b10)
        y = [1;0];
    end

end

function [uncoded_bits] = polar_SC_decoder4shaping(p, uncoded_bits, type_flag, begin_layers, end_layers)
N = length(type_flag);

soft_info = zeros(2*N-1, 1);
hard_info = zeros(2*N-1, 2);
soft_info(N:end) = log(p./(1-p));

minsum_f = @(LLR1, LLR2) sign(LLR1)*sign(LLR2)*min(abs(LLR1), abs(LLR2));
minsum_g = @(LLR1, LLR2, Bit) (1-2*Bit)*LLR1+LLR2;

for bit_index = 1:N
    %%% Calculate Soft Information Forword
    begin_layer = begin_layers(bit_index);
    for layer = begin_layer-1:-1:0
        if bit_index>1&&begin_layer-1==layer
            for index = 0:2^layer-1
                LLR1 = soft_info(2^(layer+1)+index);
                LLR2 = soft_info(2^(layer+1)+index+2^layer);
                Bit = hard_info(2^layer+index, 1);
                soft_info(2^layer+index) = minsum_g(LLR1, LLR2, Bit);
            end
        else
            for index = 0:2^layer-1
                LLR1 = soft_info(2^(layer+1)+index);
                LLR2 = soft_info(2^(layer+1)+index+2^layer);
                soft_info(2^layer+index) = minsum_f(LLR1, LLR2);
            end
        end
    end
    %%% Hard Decision By Soft Information
    if 2==type_flag(bit_index)
        hard_info(1, 1+mod(bit_index+1, 2)) = uncoded_bits(bit_index);
    else
        hard_info(1, 1+mod(bit_index+1, 2)) = soft_info(1)<0;
        uncoded_bits(bit_index) = soft_info(1)<0;
    end
    %%% Calculate Hard Information Backword
    end_layer = end_layers(bit_index);
    for layer = 0:end_layer-1
        temp_flag = (end_layer-1~=layer);
        for index = 0:2^layer-1
            hard_info(2^(layer+1)+index, 1+temp_flag) = mod(hard_info(2^layer+index, 1)+hard_info(2^layer+index, 2), 2);
            hard_info(2^(layer+1)+2^layer+index, 1+temp_flag) = hard_info(2^layer+index, 2);
        end
    end
end
end


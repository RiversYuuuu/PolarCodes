function [codeword_hat, info_hat] = polar_SC_decoder(LLR, frozen_flag, begin_layers, end_layers)
N = length(LLR);
K = sum(frozen_flag==0);
info_hat = zeros([K, 1]);
info_index = 1;

soft_info = zeros(2*N-1, 1);
hard_info = zeros(2*N-1, 2);
soft_info(N:end) = LLR;

minsum_f = @(LLR1, LLR2) sign(LLR1)*sign(LLR2)*min(abs(LLR1), abs(LLR2));
minsum_g = @(LLR1, LLR2, Bit) (1-2*Bit)*LLR1+LLR2;

for bit_index = 0:N-1
    %%% Calculate Soft Information Forword
    begin_layer = begin_layers(bit_index+1);
    for layer = begin_layer-1:-1:0
        if 0~=bit_index && begin_layer-1==layer
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
    if 1 == frozen_flag(bit_index+1)
        hard_info(1, 1+mod(bit_index, 2)) = 0;
    else
        hard_info(1, 1+mod(bit_index, 2)) = soft_info(1)<0;
        info_hat(info_index,1) =  soft_info(1)<0;
        info_index = info_index+1;
    end
    %%% Calculate Hard Information Backword
    end_layer = end_layers(bit_index+1);
    for layer = 0:end_layer-1
        temp_flag = (end_layer-1~=layer);
        for index = 0:2^layer-1
            hard_info(2^(layer+1)+index, 1+temp_flag) = mod(hard_info(2^layer+index, 1) + hard_info(2^layer+index, 2), 2);
            hard_info(2^(layer+1)+index+2^layer, 1+temp_flag) = hard_info(2^layer+index, 2);
        end
    end
end
codeword_hat = hard_info(N:end,1);
end


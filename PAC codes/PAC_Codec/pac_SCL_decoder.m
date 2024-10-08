function [codeword_hat, info_hat] = pac_SCL_decoder(LLR, conv_vec, list_size, frozen_flag, begin_layers, end_layers)
N = length(LLR);
K = sum(frozen_flag==0);
conv_len = length(conv_vec);

cache = zeros(list_size, conv_len);
path_metric = zeros(list_size, 1);
active_flag = zeros(list_size, 1);
soft_info = zeros(2*N-1, list_size);
hard_info = zeros(2*N-1, 2*list_size);
info_hat = zeros([K, list_size]);

soft_info(N:end, 1) = LLR;
active_flag(1) = 1;
info_index = 1;

minsum_f = @(LLR1, LLR2) sign(LLR1)*sign(LLR2)*min(abs(LLR1), abs(LLR2));
minsum_g = @(LLR1, LLR2, Bit) (1-2*Bit)*LLR1+LLR2;

for bit_index = 0:N-1
    %%% Calculate Soft Information Forword
    begin_layer = begin_layers(bit_index+1);
    for list_index = 1:list_size
        if active_flag(list_index) == 0
            continue
        end
        for layer = begin_layer-1:-1:0
            if 0~=bit_index && begin_layer-1==layer
                for index = 0:2^layer-1
                    LLR1 = soft_info(2^(layer+1)+index, list_index);
                    LLR2 = soft_info(2^(layer+1)+index+2^layer, list_index);
                    Bit = hard_info(2^layer+index, 2*list_index-1);
                    soft_info(2^layer+index, list_index) = minsum_g(LLR1, LLR2, Bit);
                end
            else
                for index = 0:2^layer-1
                    LLR1 = soft_info(2^(layer+1)+index, list_index);
                    LLR2 = soft_info(2^(layer+1)+index+2^layer, list_index);
                    soft_info(2^layer+index, list_index) = minsum_f(LLR1, LLR2);
                end
            end
        end
    end
    %%% Hard Decision By Soft Information
    if 1 == frozen_flag(bit_index+1)
        for list_index = 1:list_size
            if active_flag(list_index)==0
                continue
            end
            cache(list_index, :) = [0, cache(list_index, 1:end-1)];
            hard_info(1, 2*list_index-bitxor(mod(bit_index, 2), 1)) = mod(sum(conv_vec.*cache(list_index, :)), 2);
            if 0==hard_info(1, 2*list_index-bitxor(mod(bit_index, 2), 1))
                path_metric(list_index) = path_metric(list_index)+log(1+exp(-soft_info(1,list_index)));
            else
                path_metric(list_index) = path_metric(list_index)+log(1+exp(soft_info(1,list_index)));
            end
        end
    else
        current_list_size = sum(active_flag);
        %%% List is not filled
        if current_list_size < list_size
            for list_index = 1:current_list_size
                new_index = list_index+current_list_size;
                active_flag(new_index) = 1;
                soft_info(:, new_index) = soft_info(:, list_index);
                info_hat(:, new_index) = info_hat(:, list_index);
                hard_info(:, [2*new_index-1, 2*new_index]) = hard_info(:, [2*list_index-1, 2*list_index]);
                cache(list_index, :) = [0, cache(list_index, 1:end-1)];
                cache(new_index, :) = [1, cache(list_index, 1:end-1)];
                info_hat(info_index, [list_index, new_index]) = [0, 1];
                hard_info(1, [2*list_index-bitxor(mod(bit_index, 2), 1), 2*new_index-bitxor(mod(bit_index, 2), 1)]) = [mod(sum(conv_vec.*cache(list_index, :)), 2), mod(sum(conv_vec.*cache(new_index, :)), 2)];
                if 0==hard_info(1, 2*list_index-bitxor(mod(bit_index, 2), 1))
                    path_metric([list_index, new_index]) = path_metric(list_index) + [log(1 + exp(-soft_info(1,list_index))); log(1 + exp(soft_info(1,new_index)))];
                else
                    path_metric([list_index, new_index]) = path_metric(list_index) + [log(1 + exp(soft_info(1,list_index))); log(1 + exp(-soft_info(1,new_index)))];
                end
            end
        %%% List is filled
        else
            %%% Candidate Paths Increased From list_size To 2×list_size
            path_metric_double = zeros(2, list_size);
            for list_index = 1:list_size
                path_metric_double(1, list_index) = path_metric(list_index)+log(1+exp(-soft_info(1, list_index)));
                path_metric_double(2, list_index) = path_metric(list_index)+log(1+exp(soft_info(1, list_index)));
            end
            %%% Find Out the largest list_size Number of Candidate Paths
            sorted_path_metric_double = sort(path_metric_double(:));
            path_metric_middle = sorted_path_metric_double(list_size);
            equal_middle_num = sum(sorted_path_metric_double(1:list_size)==path_metric_middle);
            equal_middle_count = 0;
            compare_result = zeros(2, list_size);
            unactive_index = [];
            for list_index = 1:list_size
                for bit_value = 0:1
                    if path_metric_double(1+bit_value, list_index)<path_metric_middle
                        compare_result(1+bit_value, list_index) = 1;
                    elseif path_metric_double(1+bit_value, list_index)==path_metric_middle && equal_middle_count<equal_middle_num
                        compare_result(1+bit_value, list_index) = 1;
                        equal_middle_count = equal_middle_count+1;
                    else
                        compare_result(1+bit_value, list_index) = 0;
                    end
                end
                if compare_result(1, list_index)+compare_result(2, list_index)==0
                    unactive_index(end+1) = list_index;
                end
            end
            %%% Retain list_size Number of Candidate Paths
            for list_index = 1:list_size
                if compare_result(1, list_index)==0 && compare_result(2, list_index)==1
                    hard_info(1, 2*list_index-bitxor(mod(bit_index, 2), 1)) = 1;
                    if 1==mod([0, cache(list_index, 1:end-1)].*conv_vec, 2)
                        info_hat(info_index, list_index) = 0;
                        cache(list_index, :) = [0, cache(list_index, 1:end-1)];
                    else
                        info_hat(info_index, list_index) = 1;
                        cache(list_index, :) = [1, cache(list_index, 1:end-1)];
                    end
                    path_metric(list_index) = path_metric_double(2, list_index);
                elseif compare_result(1, list_index)==1 && compare_result(2, list_index)==0
                    hard_info(1, 2*list_index-bitxor(mod(bit_index, 2), 1)) = 0;
                    if 0==mod([0, cache(list_index, 1:end-1)].*conv_vec, 2)
                        info_hat(info_index, list_index) = 0;
                        cache(list_index, :) = [0, cache(list_index, 1:end-1)];
                    else
                        info_hat(info_index, list_index) = 1;
                        cache(list_index, :) = [1, cache(list_index, 1:end-1)];
                    end
                    path_metric(list_index) = path_metric_double(1, list_index);
                elseif compare_result(1, list_index)==1 && compare_result(2, list_index)==1
                    new_index = unactive_index(end);
                    unactive_index(end) = [];
                    soft_info(:, new_index) = soft_info(:, list_index);
                    info_hat(:, new_index) = info_hat(:, list_index);
                    cache(new_index, :) = cache(list_index, :);
                    hard_info(:, [2*new_index-1, 2*new_index]) = hard_info(:, [2*list_index-1, 2*list_index]);
                    hard_info(1, [2*list_index-bitxor(mod(bit_index, 2), 1), 2*new_index-bitxor(mod(bit_index, 2), 1)]) = [0, 1];
                    if 0==mod([0, cache(list_index, 1:end-1)].*conv_vec, 2)
                        info_hat(info_index, list_index) = 0;
                        info_hat(info_index, new_index) = 1;
                        cache(list_index, :) = [0, cache(list_index, 1:end-1)];
                        cache(new_index, :) = [1, cache(new_index, 1:end-1)];
                    else
                        info_hat(info_index, list_index) = 1;
                        info_hat(info_index, new_index) = 0;
                        cache(list_index, :) = [1, cache(list_index, 1:end-1)];
                        cache(new_index, :) = [0, cache(new_index, 1:end-1)];
                    end
                    path_metric([list_index, new_index]) = [path_metric_double(1, list_index), path_metric_double(2, list_index)];
                end
            end
        end
        info_index = info_index + 1;
    end
    %%% Calculate Hard Information Backword
    end_layer = end_layers(bit_index+1);
    for list_index = 1:list_size
        if active_flag(list_index) == 0
            continue;
        end
        for layer = 0:end_layer-1
            end_layer_flag = (end_layer-1~=layer);
            for index = 0:2^layer-1
                hard_info(2^(layer+1)+index, 2*list_index-bitxor(end_layer_flag, 1)) = mod(hard_info(2^layer+index, 2*list_index-1) + hard_info(2^layer+index, 2*list_index), 2);
                hard_info(2^(layer+1)+index+2^layer, 2*list_index-bitxor(end_layer_flag, 1)) = hard_info(2^layer+index, 2*list_index);
            end
        end
    end
end
[~, sorted_index] = sort(path_metric);
codeword_hat = hard_info(N:end,2*sorted_index-1);
info_hat = info_hat(:, sorted_index);

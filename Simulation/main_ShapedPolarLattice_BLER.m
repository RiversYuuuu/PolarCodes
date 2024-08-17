clear
clc
%%% Parameter
N = 128;
R = [0.11 0.77 0.99 1 1];
K = round(N*R);
send_power = 9;

list_size = 32;
SNR_range = 16:0.5:18.5;
max_iteration = 1e6;
max_error = 200;
patch_size = 100;
level_num = length(R);

%%% Shaping Parameter
send_prob = normpdf(0, -2^(level_num-1)+1:2^(level_num-1), sqrt(send_power));
send_prob = send_prob / sum(send_prob);
zero_prob_levels = cell(1, level_num);
prob_levels = cell(1, level_num);
for level = 1:level_num
    flag = mod(1-2^(level_num-1):2^(level_num-1),2^level);
    zero_prob = zeros(1, 2^(level-1));
    prob = zeros(1, 2^(level-1));
    for i = 0:2^(level-1)-1
        prob(i+1) = (sum(send_prob(flag==i+2^(level-1)))+sum(send_prob(flag==i)));
        zero_prob(i+1) = sum(send_prob(flag==i))/(sum(send_prob(flag==i+2^(level-1)))+sum(send_prob(flag==i)));
    end
    zero_prob_levels{1, level} = zero_prob;
    prob_levels{1, level} = prob;
end

entropy = zeros(1, level_num);
for level = 1:level_num
    entropy(level) = sum(get_entropy(zero_prob_levels{1, level}).*prob_levels{1, level});
end
K_shaping = round(N*(1-entropy));
K_info = round(N*R)-K_shaping;
K_frozen = N-K_info-K_shaping;

%%% Construction For Transmission
construct_SNR = SNR_range(end);
construct_sigma = sqrt(send_power/10^(construct_SNR/10));
frozen_flag = zeros(level_num, N);
sorted_indice4frozen = zeros(N, level_num);
error_bound4frozen = zeros(N, level_num);
for level = 1 : level_num
    [error_bound4frozen(:, level), sorted_indice4frozen(:, level)] = polar_degrademerge_ModAWGN(N, construct_sigma, level, ones(1, 2^level_num)/2^level_num);
    frozen_flag(level, sorted_indice4frozen(end+1-K_frozen(level):end, level)) = 1;
end

%%% Construction For Shaping
shaping_flag = zeros(level_num, N);
sorted_indice4shaping = zeros(N, level_num);
error_bound4shaping = zeros(N, level_num);
for level = 1:level_num
    Pary = [prob_levels{1, level}; zero_prob_levels{1, level}];
    [~, sorted_indice4shaping(:, level)] = polar_degrademerge_BSC(N, Pary, 16);
    shaping_flag(level, sorted_indice4shaping(1:K_shaping(level), level)) = 1;
end

%%% Simulation
[begin_layers, end_layers] = polar_decode_prepare(N);
type_flag = frozen_flag+2*shaping_flag;
error_count = zeros(level_num+1, length(SNR_range));
BLER = zeros(size(SNR_range));
lambda_sum = [];
for SNR_index = 1:length(SNR_range)
    SNR = SNR_range(SNR_index);
    noise_sigma = sqrt(send_power/10^(SNR/10));
    msg_length = [0, 0];
    for patch_index = 1:max_iteration/patch_size
        info_hat = cell(1, patch_size);
        lambda_hat = zeros(N, patch_size);
        info_patch = cell(1, patch_size);
        lambda = zeros(N, patch_size);
        parfor iter = 1:patch_size
            %%% Encoding
            info = cell(1, level_num);
            for level = 1 : level_num
                info{1, level} = rand(K_info(level), 1)>0.5;
            end
            info_patch{iter} = info;
            lambda(:, iter) = polar_lattice_encoder4shaping(info, type_flag, zero_prob_levels, begin_layers, end_layers);
            %%% Transmission
            recv_signal = lambda(:, iter)+noise_sigma*randn(N, 1);
            %%% Decoding
            [lambda_hat(:, iter), info_hat{iter}] = polar_lattice_decoder4shaping(recv_signal, noise_sigma, zero_prob_levels, type_flag, list_size, begin_layers, end_layers);
        end
        %%% Statistics
        for iter = 1:patch_size
            for level = 1:level_num
                error_count(level, SNR_index) = error_count(level, SNR_index) + ~isequal(info_patch{iter}{level}, info_hat{iter}{level});
            end
            error_count(end, SNR_index) = error_count(end, SNR_index) + ~isequal(info_patch{iter}, info_hat{iter});
        end
        fprintf(repmat('\b', 1, 1+sum(msg_length)));
        msg1 = ['SNR = ',num2str(SNR,'%.2f'),' dB, Run ', num2str(patch_index*patch_size    ), '/', num2str(max_iteration), ' times'];
        msg2 = ['Error ', num2str(error_count(end, SNR_index)), ' times, ', 'BLER = ', num2str(error_count(end, SNR_index)/iter/patch_index, '%.2e')];
        msg_length = [length(msg1), length(msg2)];
        fprintf([msg1, '\n', msg2]);
        if error_count(end, SNR_index)>=max_error||iter==max_iteration
            fprintf('\n\n\n');
            break;
        end
    end
    BLER(SNR_index) = error_count(SNR_index)/iter/patch_size;
end



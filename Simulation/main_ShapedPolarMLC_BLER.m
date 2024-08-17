clear
clc
%%% Parameter
N = 4096;
R = [0.04 0.31 0.70 0.96 0.99 1];
K = round(N*R);
modulate_type = '64Spiral';

list_size = 1;
SNR_range = 13.5:0.25:14.25;
max_iteration = 1e6;
max_error = 200;
patch_size = 100;

%%% Modulation
[send_set, send_prob_set] = get_send_set(modulate_type, 1);
level_num = log2(length(send_set));
send_power = sum(abs(send_set).^2.*send_prob_set);
zero_prob_levels = cell(1, level_num);
prob_levels = cell(1, level_num);
for level=1:level_num
    step = 2^level;
    for i=1:2^(level-1)
        prob0 = sum(send_prob_set(i:step:end));
        prob1 = sum(send_prob_set(i+2^(level-1):step:end));
        zero_prob_levels{level}(i) = prob0/(prob0+prob1);
        prob_levels{level}(i) = prob0+prob1;
    end
end

%%% Shaping Parameter
zero_prob_levels = cell(1, level_num);
prob_levels = cell(1, level_num);
for level = 1:level_num
    zero_prob = zeros(1, 2^(level-1));
    prob = zeros(1, 2^(level-1));
    step = 2^level;
    for i = 0:2^(level-1)-1
        zero_index = find(mod(0:2^level_num-1, 2^level)==i);
        one_index = find(mod(0:2^level_num-1, 2^level)==i+2^(level-1));
        prob(i+1) = sum(send_prob_set(one_index)) + sum(send_prob_set(zero_index));
        zero_prob(i+1) = sum(send_prob_set(zero_index))/prob(i+1);
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
construct_sigma = sqrt(send_power/(2*10^(construct_SNR/10)));
frozen_flag = zeros(level_num, N);
[error_bound4frozen, sorted_indice4frozen] = polar_degrademerge_MLC(N, construct_sigma, send_set, send_prob_set);
for level = 1:level_num
    frozen_flag(level, sorted_indice4frozen(K(level)+1:end, level)) = 1;
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
for SNR_index = 1:length(SNR_range)
    SNR = SNR_range(SNR_index);
    noise_sigma = sqrt(send_power/(2*10^(SNR/10)));
    msg_length = [0, 0];
    for patch_index = 1:max_iteration/patch_size
        info_hat = cell(1, patch_size);
        parfor iter = 1:patch_size
            %%% Encoding
            info = cell(1, level_num);
            for level = 1 : level_num
                info{1, level} = rand(K(level), 1)>0.5;
            end
            send_signal = polar_MLC_encoder4shaping(info, send_set, frozen_flag);
            %%% Transmission
            recv_signal = send_signal+noise_sigma*(randn(N, 1)+randn(N, 1)*1i);
            %%% Decoding
            [~, info_hat{iter}] = polar_MLC_decoder(send_set, send_prob_set, recv_signal, noise_sigma, frozen_flag, list_size, begin_layers, end_layers);
        end
        %%% Statistics
        for iter = 1:patch_size
            for level = 1:level_num
                error_count(level, SNR_index) = error_count(level, SNR_index) + ~isequal(info_patch{iter}{level}, info_hat{iter}{level});
            end
            error_count(end, SNR_index) = error_count(end, SNR_index) + ~isequal(info_patch{iter}, info_hat{iter});
        end
        fprintf(repmat('\b', 1, 1+sum(msg_length)));
        msg1 = ['SNR = ',num2str(SNR, '%.2f'),' dB, Run ', num2str(patch_size*patch_index), '/', num2str(max_iteration), ' times'];
        msg2 = ['Error ', num2str(error_count(end, SNR_index)), ' times, ', 'BLER = ', num2str(error_count(end, SNR_index)/iter/patch_index, '%.2e')];
        msg_length = [length(msg1), length(msg2)];
        fprintf([msg1, '\n', msg2]);
        if error_count(end, SNR_index)>=max_error||patch_index==max_iteration/patch_size
            fprintf('\n\n\n');
            break;
        end
    end
    BLER(SNR_index) = error_count(SNR_index)/iter/patch_size;
end
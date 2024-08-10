clear
clc
%%% Parameter
N = 1024;
R = [0.42 0.95 1 1 1];
K = round(N*R);
send_power = 9;

list_size = 1;
SNR_range = 11:0.5:15;
max_iteration = 1e6;
max_error = 200;
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
frozen_flag = zeros(N, level_num);
sorted_indice4frozen = zeros(N, level_num);
error_bound4frozen = zeros(N, level_num);
for level = 1 : level_num
    [error_bound4frozen(:, level), sorted_indice4frozen(:, level)] = polar_degrademerge_ModAWGN(N, construct_sigma, level, ones(1, 2^level_num)/2^level_num);
    frozen_flag(sorted_indice4frozen(end+1-K_frozen(level):end, level), level) = 1;
end

%%% Construction For Shaping
shaping_flag = zeros(N, level_num);
sorted_indice4shaping = zeros(N, level_num);
error_bound4shaping = zeros(N, level_num);
for level = 1:level_num
    Pary = [prob_levels{1, level}; zero_prob_levels{1, level}];
    [~, sorted_indice4shaping(level, :)] = polar_degrademerge_BSC(N, Pary, 16);
    shaping_flag(sorted_indice4shaping(1:K_shaping(level), level), level) = 1;
end

%%% Simulation
[begin_layers, end_layers] = polar_decode_prepare(N);
type_flag = frozen_flag+2*shaping_flag;
error_count = zeros(level_num+1, length(SNR_range));
BLER = zeros(size(SNR_range));
for SNR_index = 1:length(SNR_range)
    SNR = SNR_range(SNR_index);
    noise_sigma = sqrt(send_power/10^(SNR/10));
    msg_length = [0, 0];
    for iter = 1:max_iteration
        %%% Encoding
        info = cell(1, level_num);
        for level = 1 : level_num
            info{1, level} = rand(K_info(level), 1)>0.5;
        end
        lambda = polar_lattice_encoder4shaping(info, type_flag, zero_prob_levels, begin_layers, end_layers);
        %%% Transmission
        recv_signal = lambda+noise_sigma*randn(N, 1);
        %%% Decoding
        info_hat = cell(1, level_num);
        LLR = zeros(N, level_num);
        for level = 1:level_num
            p0 = zeros(N, 1);
            p1 = zeros(N, 1);
            mod_recv_signal = mod(recv_signal, 2);
            sigma_level = noise_sigma/2^(level-1);
            for bias = -8:2:8
                p0 = p0+normpdf(mod_recv_signal, 0+bias, sigma_level);
                p1 = p1+normpdf(mod_recv_signal, 1+bias, sigma_level);
            end
            LLR(:, level) = log(p0./p1);
            [~, info_hat{level}] = polar_SCL_decoder(LLR(:, level), list_size, frozen_flag(:, level), begin_layers, end_layers);
            recv_signal = 0.5*(recv_signal-polar_encoder4lattice(info_hat{level}, frozen_flag(:, level)));
            error_count(level, SNR_index) = error_count(level, SNR_index)+~isequal(info{level}, info_hat{level});
        end
        %%% Statistics
        error_count(end, SNR_index) = error_count(end, SNR_index) + ~isequal(info, info_hat);
        if mod(iter, 100)==0
            fprintf(repmat('\b', 1, 1+sum(msg_length)));
            msg1 = ['SNR = ',num2str(SNR,'%.2f'),' dB, Run ', num2str(iter), '/', num2str(max_iteration), ' times'];
            msg2 = ['Error ', num2str(error_count(end, SNR_index)), ' times, ', 'BLER = ', num2str(error_count(end, SNR_index)/iter, '%.2e')];
            msg_length = [length(msg1), length(msg2)];
            fprintf([msg1, '\n', msg2]);
            if error_count(end, SNR_index)>=max_error||iter==max_iteration
                fprintf('\n\n\n');
                break;
            end
        end
    end
    BLER(SNR_index) = error_count(SNR_index)/iter;
end



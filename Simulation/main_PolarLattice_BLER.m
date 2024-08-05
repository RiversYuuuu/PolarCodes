clear
clc
%%% Parameter
N = 1024;
R = [0.22, 0.91];
K = round(N*R);

list_size = 1;
VNR_range = 0:0.5:2.5;
max_iteration = 1e6;
max_error = 200;
V_lamda = 2^(N*2-N*sum(R));
level_num = length(R);
%%% Construction
construct_VNR = VNR_range(end);
construct_sigma = 1/sqrt(10^(construct_VNR/10)*(2*pi*exp(1)))*(V_lamda^(1/N));
frozen_flag = ones(N, level_num);
sorted_indice = zeros(N, level_num);
error_bound = zeros(N, level_num);
for level = 1 : level_num
    [error_bound(:, level), sorted_indice(:, level)] = polar_degrademerge_ModAWGN(N, construct_sigma, level, ones(1, 2^level_num)/2^level_num);
    frozen_flag(sorted_indice(1:K(level), level), level) = 0;
end

[begin_layers, end_layers] = polar_decode_prepare(N);
error_count = zeros(level_num+1, length(VNR_range));
BLER = zeros(size(VNR_range));
for VNR_index = 1:length(VNR_range)
    VNR = VNR_range(VNR_index);
    noise_sigma = 1 / sqrt(10^(VNR/10)*(2*pi*exp(1))) * (V_lamda^(1/N));
    msg_length = [0, 0];
    for iter = 1:max_iteration
        %%% Encoding
        info = cell(1, level_num);
        for level = 1 : level_num
            info{1, level} = rand(K(level), 1)>0.5;
        end
        lambda = polar_lattice_encoder(info, frozen_flag);
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
            error_count(level, VNR_index) = error_count(level, VNR_index)+~isequal(info{level}, info_hat{level});
        end
        %%% Statistics
        error_count(end, VNR_index) = error_count(end, VNR_index) + ~isequal(info, info_hat);
        if mod(iter, 100)==0
            fprintf(repmat('\b', 1, 1+sum(msg_length)));
            msg1 = ['VNR = ',num2str(VNR,'%.2f'),' dB, Run ', num2str(iter), '/', num2str(max_iteration), ' times'];
            msg2 = ['Error ', num2str(error_count(end, VNR_index)), ' times, ', 'BLER = ', num2str(error_count(end, VNR_index)/iter, '%.2e')];
            msg_length = [length(msg1), length(msg2)];
            fprintf([msg1, '\n', msg2]);
            if error_count(end, VNR_index)>=max_error||iter==max_iteration
                fprintf('\n\n\n');
                break;
            end
        end
    end
    BLER(VNR_index) = error_count(VNR_index)/iter;
end



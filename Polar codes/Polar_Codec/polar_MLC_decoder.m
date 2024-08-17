function [codeword_hat, info_hat] = polar_MLC_decoder(send_set, send_prob_set, recv_signal, noise_sigma, frozen_flag, list_size, begin_layers, end_layers)
[level_num, N] = size(frozen_flag);
LLR = zeros(N, level_num);
dim = 2-(sum(abs(imag(send_set)))==0);
info_hat = cell(1, level_num);
codeword_hat = zeros(N, level_num);
for level = 1:level_num
    step = 2^level;
    index = 1+bi2de(codeword_hat);
    for i = 1:N
        zero_index = index(i):step:length(send_set);
        one_index = index(i)+2^(level-1):step:length(send_set);
        p0 = send_prob_set(zero_index)*get_AWGN_send_prob(send_set(zero_index), recv_signal(i), noise_sigma, dim);
        p1 = send_prob_set(one_index)*get_AWGN_send_prob(send_set(one_index), recv_signal(i), noise_sigma, dim);
        LLR(i, level) = log(p0/p1);
    end
    [codeword_hat(:, level), uncoded_bits] = polar_SCL_decoder_mex(LLR(:, level), list_size, frozen_flag(level, :), begin_layers, end_layers);
    info_hat{level} = uncoded_bits(frozen_flag(level, :)==0);
end
end
function [lambda_hat, info_hat] = polar_lattice_decoder4shaping(recv_signal, noise_sigma, zero_prob_levels, type_flag, list_size, begin_layers, end_layers)
[level_num, N] = size(type_flag);
LLR = zeros(N, level_num);
info_hat = cell(1, level_num);
lambda_hat = zeros(N, 1);
for level = 1:level_num
    mod_recv_signal = mod(recv_signal, 2);
    sigma_level = noise_sigma/2^(level-1);
    bias = -2^level_num:2:2^level_num;
    zero_prob_index = 1+mod(lambda_hat, 2^(level-1));
    reverse_flag = floor(mod(lambda_hat, 2^(level))/2^(level-1));
    p = reshape(zero_prob_levels{level}(zero_prob_index), [], 1);
    p(reverse_flag==1) = 1-p(reverse_flag==1);
    p0 = p.*sum(normpdf(mod_recv_signal, 0+bias, sigma_level), 2);
    p1 = (1-p).*sum(normpdf(mod_recv_signal, 1+bias, sigma_level), 2);
    LLR(:, level) = log(p0./p1);
    [~, uncoded_bits] = polar_SCL_decoder_mex(LLR(:, level), list_size, type_flag(level, :), begin_layers, end_layers);
    info_hat{level} = uncoded_bits(type_flag(level, :)==0);
    codeword = polar_encoder4lattice(uncoded_bits);
    recv_signal = 0.5*(recv_signal-codeword);
    lambda_hat = lambda_hat+2^(level-1)*codeword;
end
end
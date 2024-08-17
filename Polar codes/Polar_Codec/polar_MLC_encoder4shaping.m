function [send_signal] = polar_MLC_encoder4shaping(info, send_set, type_flag, zero_prob_levels, begin_layers, end_layers)
[level_num, N] = size(frozen_flag);
codeword = zeros(N, level_num);
for level = 1:level_num
    uncoded_bits = zeros(N, 1);
    uncoded_bits(type_flag(level, :)==0) = info{1, level};
    zero_prob_index = bi2de(codeword);
    p = zero_prob_levels{level}(zero_prob_index);
    [~, codeword(:, level)] = polar_SC_decoder4shaping(p, uncoded_bits, type_flag(level, :), begin_layers, end_layers);
end
send_signal = reshape(send_set(1+bi2de(codeword)), [], 1);
end
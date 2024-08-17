function [lambda] = polar_lattice_encoder4shaping(info, type_flag, zero_prob_levels, begin_layers, end_layers)
[level_num, N] = size(type_flag);
lambda = zeros(N, 1);
for level = 1:level_num
    uncoded_bits = zeros(N, 1);
    uncoded_bits(type_flag(level, :)==0) = info{1, level};
    zero_prob_index = 1+mod(lambda, 2^(level-1));
    reverse_flag = floor(mod(lambda, 2^(level))/2^(level-1));
    p = zero_prob_levels{level}(zero_prob_index);
    p(reverse_flag==1) = 1-p(reverse_flag==1);
    uncoded_bits = polar_SC_decoder4shaping(p, uncoded_bits, type_flag(level, :), begin_layers, end_layers);
    lambda = lambda+2^(level-1)*polar_encoder4lattice(uncoded_bits);
end
lambda = mod(lambda, 2^level_num);
minus_index = find(lambda>2^(level_num-1));
lambda(minus_index) = lambda(minus_index) - 2^level_num;
end
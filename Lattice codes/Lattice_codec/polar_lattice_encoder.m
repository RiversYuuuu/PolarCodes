function [lambda] = polar_lattice_encoder(info, frozen_flag)
[N, level_num] = size(frozen_flag);
lambda = zeros(N, 1);
for level = 1:level_num
    uncoded_bits = zeros(N, 1);
    uncoded_bits(frozen_flag==0) = info{1, level};
    lambda = lambda+2^(level-1)*polar_encoder4lattice(uncoded_bits);
end
lambda = mod(lambda, 2^level_num);
minus_index = find(lambda>2^(level_num-1));
lambda(minus_index) = lambda(minus_index) - 2^level_num;
end
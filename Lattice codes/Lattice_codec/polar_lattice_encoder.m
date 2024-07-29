function [lambda] = polar_lattice_encoder(info, frozen_flag)
[N, level_num] = size(frozen_flag);
lambda = zeros(N, 1);
for level = 1:level_num
    lambda = lambda+2^(level-1)*polar_encoder4lattice(info{1, level}, frozen_flag(:, level));
end
lambda = mod(lambda, 2^level_num);
minus_index = find(lambda>2^(level_num-1));
lambda(minus_index) = lambda(minus_index) - 2^level_num;
end
function [error_bound, sorted_indices] = polar_degrademerge_ModAWGN(N, noise_sigma, level, send_prob)
quan_level = 32;
channel_type_num = 2^(level-1);
noise_sigma = noise_sigma/channel_type_num;
level_num = log2(length(send_prob));
%%% Channel Types
flag = mod(1-2^(level_num-1):2^(level_num-1), 2^level);
back_level_prob = zeros(1, channel_type_num);
back_level_prob0 = zeros(1, channel_type_num);
for i = 0:channel_type_num-1
    prob0 = sum(send_prob(flag==i));
    prob1 = sum(send_prob(flag==i+channel_type_num));
    back_level_prob(i+1) = prob0+prob1;
    back_level_prob0(i+1) = prob0/back_level_prob(i+1);
end
%%% Channel Quantization
p0 = zeros(channel_type_num, quan_level);
p1 = zeros(channel_type_num, quan_level);
y = linspace(0, 2, quan_level);
for i = 1:channel_type_num
    for bias = -8:2:8
        p0(i, :) = p0(i, :)+back_level_prob(i)*back_level_prob0(i)*normpdf(y+bias, 0, noise_sigma);
        p1(i, :) = p1(i, :)+back_level_prob(i)*(1-back_level_prob0(i))*normpdf(y+bias, 1, noise_sigma);
    end
end
p0 = reshape(p0, 1, []);
p1 = reshape(p1, 1, []);
P_array = [(p0+p1)/sum(p0+p1); min([p0; p1])./(p0+p1)];
[error_bound, sorted_indices] = polar_degrademerge_BSC_Z(N, P_array, quan_level);
end
function [error_bound, sorted_indices] = polar_degrademerge_MLC(N, construct_sigma, send_set, send_prob_set)
%%% Parameter
quan_level = 16;
level_num = log2(length(send_set));
error_bound = zeros(N, level_num);
sorted_indices = zeros(N, level_num);

%%% Channel Quantization
dim = 2-(sum(abs(imag(send_set)))==0);
quan_max = max(real(send_set))+3*construct_sigma;
quan_range = linspace(quan_max, -quan_max, quan_level);
recv_signal = zeros(1, length(quan_range)^dim);
recv_index = 1;
for recv_real = quan_range
    recv_signal(recv_index) = recv_real;
    for recv_imag = quan_range
        if dim==1
            recv_index = recv_index + 1;
            break
        end
        recv_signal(recv_index) = recv_real + recv_imag * 1i;
        recv_index = recv_index + 1;
    end
end
for level = 1:level_num
    step = 2^level;
    p0 = [];
    p1 = [];
    for i = 1:2^(level-1)
        zero_send_signal = send_set(i:step:end);
        one_send_signal = send_set(i+2^(level-1):step:end);
        p0 = [p0, send_prob_set(i:step:end)*get_AWGN_send_prob(zero_send_signal, recv_signal, construct_sigma, dim)];
        p1 = [p1, send_prob_set(i+2^(level-1):step:end)*get_AWGN_send_prob(one_send_signal, recv_signal, construct_sigma, dim)];
    end
    P_array = [(p0+p1)/sum(p0+p1); min([p0; p1])./(p0+p1)];
    [error_bound(:, level), sorted_indices(:, level)] = polar_degrademerge_BSC(N, P_array, quan_level);
    fprintf("Construction of %d-th level is done.\n",level);
end
end
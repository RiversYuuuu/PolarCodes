function [best_order] = get_maxVar_order(send_set, send_prob_set, SNR)
M = length(send_set);
order_set = perms(1:log2(M));
set_num = size(order_set, 1);
mutual_info_var = zeros(1, set_num);
parfor order_index = 1:set_num
    index = de2bi(0:M-1);
    index = bi2de(index(:, order_set(order_index, :)));
    [~, mutual_info_var(1, order_index), ~] = get_mutual_info(send_set(index+1), send_prob_set(index+1), SNR);
end
[~, sorted_index] = sort(mutual_info_var);
best_order = order_set(sorted_index(end), :);
end
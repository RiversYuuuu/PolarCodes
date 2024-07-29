function [error_bound, sorted_indices] = polar_degrademerge_BSC_Z(N, P_array, quan_level)
n = log2(N);
[P_array] = degrade_merge(P_array, quan_level);
Z0 = 2*sum(P_array(1,:).*sqrt(P_array(2,:).*(1-P_array(2,:))));
Z = zeros(N, 1);
error_bound = zeros(1, N);
for bit_index = 1:N
    binary_array = de2bi(bit_index-1,log2(N), 'left-msb');
    Z(bit_index) = Z0;
    Q_array = P_array;
    for binart_index = 1:n
        if binary_array(binart_index) == 0
            [~,W] = WPlusMinus(Q_array);
            Z_W = 2*sum(W(1,:).*sqrt(W(2,:).*(1-W(2,:))));
            Z(bit_index) = min([2*Z(bit_index)-Z(bit_index)^2, Z_W]);
        else
            [W,~] = WPlusMinus(Q_array);
            Z(bit_index) = Z(bit_index)^2;
        end
        Q_array = degrade_merge(W, quan_level);
    end
    error_bound(bit_index) = min([sum(Q_array(1,:).*Q_array(2,:)), Z(bit_index)]);
end
[~, sorted_indices] = sort(error_bound);











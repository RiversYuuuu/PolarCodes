function [error_prob_bound, sorted_indices] = polar_degrademerge_construction(N, P_array)
quantization_level = 16;
n = log2(N);
[P_array] = P_merge(P_array, quantization_level);
channel_layer = cell(1, 1);
channel_layer{1} = P_array;
error_prob_bound = zeros(1, N);
for layer = 1:n
    channel_next_layer = cell(1,2^layer);
    for index = 1:2^(layer-1)
        channel = channel_layer{index};
        [P_up, P_down] = channel_polarization(channel);
        P_up = P_merge(P_up, quantization_level);
        P_down = P_merge(P_down, quantization_level);
        Pe_up = sum(P_up(1,:).*P_up(2,:));
        Pe_down = sum(P_down(1,:).*P_down(2,:));
        channel_next_layer{index} = P_down;
        channel_next_layer{2^(layer-1)+index} = P_up;
        if layer==n
            error_prob_bound([index, 2^(layer-1)+index]) = [Pe_down, Pe_up];
        end
    end
    channel_layer = channel_next_layer;
end
[~, sorted_indices] = sort(error_prob_bound);











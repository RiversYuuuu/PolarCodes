function [sorted_indices] = pac_WeightedSum_construction(conv_vec, error_bound)
N = length(error_bound);
n = log2(N);
conv_len = length(conv_vec);

info_flag = zeros(1,N);
sorted_indices = zeros(1, N);

bit_weight = 2.^sum(de2bi(0:N-1),2);
current_index = 1;

%%% Estimate Sub-channels Reliability
omega_i = ceil((1-log2(1+error_bound)) / 0.1);
for weight = n:-1:0
    if nchoosek(n, weight) == 1
        selected_indices = find(bit_weight == 2^weight);
        sorted_indices(current_index) = selected_indices;
        info_flag(selected_indices) = 1;
        current_index = current_index+1;
    else
        candidate_indices = find(bit_weight == 2^weight);
        for candidate_index = 1:length(candidate_indices)
            %%% Estimate Sub-channels Utilization
            tau_i = zeros(1,N);
            for bit_index = 1:N
                for conv_index = 1:conv_len
                    if (bit_index-conv_index+1)>0
                        tau_i(bit_index) = tau_i(bit_index)+info_flag(bit_index-conv_index+1)*conv_vec(conv_index);
                    end
                end
            end
            %%% Jointly Consider Sub-channels Utilization and Reliability
            theta_i = zeros(1,N);
            for bit_index = 1:N
                for conv_index = 1:conv_len
                    if (bit_index+conv_index-1)<=N
                        theta_i(bit_index) = theta_i(bit_index)+conv_vec(conv_index)*omega_i(bit_index+conv_index-1)/(tau_i(bit_index+conv_index-1)+1);
                    end
                end
            end
            %%% Select Indices With Smallest Theta
            candidate_theta = theta_i(candidate_indices);
            max_candidate_index = find(candidate_theta == max(candidate_theta), 1);
            new_indice = candidate_indices(max_candidate_index);
            info_flag(new_indice) = 1;
            sorted_indices(1, current_index) = new_indice;
            current_index = current_index + 1;
            candidate_indices(max_candidate_index) = [];
        end
    end
end
end
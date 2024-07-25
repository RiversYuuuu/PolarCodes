function [Q_array] = degrade_merge(P_array, quantization_level)
%%% Lossless Merge
P_array(:, P_array(1,:)==0) = [];
P_unique = unique(P_array(2,:),'sorted');
unique_num = length(P_unique);
Q_array = zeros(2, unique_num);
for index = 1:unique_num
    same_flag = P_array(2,:)==P_unique(index);
    Q_array(1, index) = sum(P_array(1,same_flag));
    Q_array(2, index) = P_unique(index);
end
%%% Degrade Merge
quantization_num = length(Q_array);
Q_array_I = sqrt(Q_array(2, :).*(1-Q_array(2, :)));
while quantization_num>quantization_level
    P_sum = zeros(1, quantization_num-1);
    merged_BSC_P = zeros(1, quantization_num-1);
    merged_BSC_I = zeros(1, quantization_num-1);
    Delta_I = zeros(1, quantization_num-1);
    for index = 1:quantization_num-1
        P_sum(index) = Q_array(1, index) + Q_array(1, index+1);
        merged_BSC_P(index) = Q_array(1, index)/P_sum(index)*Q_array(2, index)+...
                       Q_array(1, index+1)/P_sum(index)*Q_array(2, index+1);
        merged_BSC_I(index) = sqrt(merged_BSC_P(index)*(1-merged_BSC_P(index)));
        Delta_I(index) = P_sum(index)*merged_BSC_I(index)-Q_array(1,index)*Q_array_I(index)-Q_array(1,index+1)*Q_array_I(index+1);
    end
    min_index = find(Delta_I == min(Delta_I));
    min_index = min_index(1);
    
    Q_array(1, min_index) = Q_array(1, min_index)+Q_array(1, min_index+1);
    Q_array(2, min_index) = merged_BSC_P(min_index);
    Q_array(:, min_index+1)=[];
    Q_array_I(min_index) = merged_BSC_I(min_index);
    Q_array_I(min_index+1) = [];
    
    quantization_num = quantization_num-1;
end

end


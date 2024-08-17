clear
clc
%%% Parameter
N = 1024;
R = [0.22, 0.91];
K = round(N*R);

list_size = 1;
VNR_range = 0:0.5:2.5;
max_iteration = 1e6;
patch_size = 100;
max_error = 200;
V_lamda = 2^(N*2-N*sum(R));
level_num = length(R);

%%% Construction
load('data\constructData_polarLattice1024_VNR2.mat')
% construct_VNR = 2;
% construct_sigma = 1/sqrt(10^(construct_VNR/10)*(2*pi*exp(1)))*(V_lamda^(1/N));
% frozen_flag = ones(level_num, N);
% sorted_indice = zeros(N, level_num);
% error_bound = zeros(N, level_num);
% for level = 1 : level_num
%     [error_bound(:, level), sorted_indice(:, level)] = polar_degrademerge_ModAWGN(N, construct_sigma, level, ones(1, 2^level_num)/2^level_num);
%     frozen_flag(level, sorted_indice(1:K(level), level)) = 0;
% end

%%% Simulation
[begin_layers, end_layers] = polar_decode_prepare(N);
error_count = zeros(level_num+1, length(VNR_range));
BLER = zeros(size(VNR_range));
for VNR_index = 1:length(VNR_range)
    VNR = VNR_range(VNR_index);
    noise_sigma = 1 / sqrt(10^(VNR/10)*(2*pi*exp(1))) * (V_lamda^(1/N));
    msg_length = [0, 0];
    for patch_index = 1:max_iteration/patch_size
        info_hat = cell(1, patch_size);
        lambda_hat = zeros(N, patch_size);
        info_patch = cell(1, patch_size);
        lambda = zeros(N, patch_size);
        parfor iter = 1:patch_size
            %%% Encoding
            info = cell(1, level_num);
            for level = 1 : level_num
                info{1, level} = rand(K(level), 1)>0.5;
            end
            info_patch{iter} = info;
            lambda(:, iter) = polar_lattice_encoder(info, frozen_flag);
            %%% Transmission
            recv_signal = lambda(:, iter)+noise_sigma*randn(N, 1);
            %%% Decoding
            [lambda_hat(:, iter), info_hat{iter}] = polar_lattice_decoder(recv_signal, noise_sigma, frozen_flag, list_size, begin_layers, end_layers);
        end
        %%% Statistics
        for iter = 1:patch_size
            for level = 1:level_num
                error_count(level, VNR_index) = error_count(level, VNR_index) + ~isequal(info_patch{iter}{level}, info_hat{iter}{level});
            end
            error_count(end, VNR_index) = error_count(end, VNR_index) + ~isequal(info_patch{iter}, info_hat{iter});
        end
        fprintf(repmat('\b', 1, 1+sum(msg_length)));
        msg1 = ['VNR = ',num2str(VNR, '%.2f'),' dB, Run ', num2str(patch_size*patch_index), '/', num2str(max_iteration), ' times'];
        msg2 = ['Error ', num2str(error_count(end, VNR_index)), ' times, ', 'BLER = ', num2str(error_count(end, VNR_index)/iter/patch_index, '%.2e')];
        msg_length = [length(msg1), length(msg2)];
        fprintf([msg1, '\n', msg2]);
        if error_count(end, VNR_index)>=max_error||patch_index==max_iteration/patch_size
            fprintf('\n\n\n');
            break;
        end
    end
    BLER(VNR_index) = error_count(VNR_index)/iter/patch_size;
end
    


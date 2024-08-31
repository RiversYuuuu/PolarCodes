clear
clc
%%% Parameter
N = 1024;
R = 0.5;
K = round(N*R);
list_size = 8;
EbN0_range = 0:0.25:2.5;
max_iteration = 1e6;
patch_size = 100;
max_error = 200;

%%% Construction
construct_EbN0 = 2;
construct_sigma = sqrt(1/R/10^(construct_EbN0/10)/2);
[error_bound, sorted_indices] = polar_degrademerge_BIAWGN(N, construct_sigma);
frozen_flag = ones(1, N);
frozen_flag(sorted_indices(1:K)) = 0;

[begin_layers, end_layers] = polar_decode_prepare(N);
error_count = zeros(size(EbN0_range));
BLER = zeros(size(EbN0_range));
for EbN0_index = 1:length(EbN0_range)
    EbN0 = EbN0_range(EbN0_index);
    noise_sigma = sqrt(1/R/10^(EbN0_range(EbN0_index)/10)/2);
    msg_length = [0, 0];
    info = rand(K, patch_size) >= 0.5;
    info_hat = zeros(K, patch_size);
    for patch_index = 1:max_iteration/patch_size
        parfor iter = 1:patch_size
            %%% Encoding
            codeword = polar_encoder_channel(info(:, iter), frozen_flag);
            %%% Transmission
            send_signals = 1-2*codeword;
            recv_signal = send_signals+noise_sigma*randn(N, 1);
            %%% Decoding
            LLR = 2*recv_signal/noise_sigma^2;
%             [~, uncoded_bits] = polar_SCL_decoder(LLR, list_size, frozen_flag, begin_layers, end_layers);
            [~, uncoded_bits] = polar_SCL_decoder_mex(LLR, list_size, frozen_flag, begin_layers, end_layers);
            info_hat(:, iter) = uncoded_bits(frozen_flag==0);
        end
        %%% Statistics
        error_count(EbN0_index) = error_count(EbN0_index)+sum(sum(info~=info_hat)~=0);
        fprintf(repmat('\b', 1, 1+sum(msg_length)));
        msg1 = ['EbN0 = ',num2str(EbN0,'%.2f'),' dB, Run ', num2str(patch_index*patch_size), '/', num2str(max_iteration), ' times'];
        msg2 = ['Error ', num2str(error_count(EbN0_index)), ' times, ', 'BLER = ', num2str(error_count(EbN0_index)/patch_index/patch_size, '%.2e')];
        msg_length = [length(msg1), length(msg2)];
        fprintf([msg1, '\n', msg2]);
        if error_count(EbN0_index)>=max_error||patch_index==max_iteration/patch_size
            fprintf('\n\n\n');
            break;
        end
    end
    BLER(EbN0_index) = error_count(EbN0_index)/patch_index/patch_size;
end

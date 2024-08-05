clear
clc
%%% Parameter
N = 2^10;
Ber_p = 0.11;
rate = 0.6;
K = round(N*rate);
max_iterasion = 1e4;


%%% Construction
high_entropy_flag = ones(N, 1);
[error_bound, sorted_indices] = polar_degrademerge_BSC_Z(N, [1; Ber_p], 16);
high_entropy_flag(sorted_indices(1:N-K)) = 0;

[begin_layers, end_layers] = polar_decode_prepare(N);
msg_length = 0;
bit_error_count = 0;
for iter=1:max_iterasion
    %%% Encoding
    info = rand(N, 1)>Ber_p;
    codeword = polar_encoder_source(info);
    high_entropy_bits = codeword(high_entropy_flag==1);
    
    %%% Decoding
    [info_hat] = polar_SC_decompressor(Ber_p, high_entropy_bits, high_entropy_flag, begin_layers, end_layers);
    
    %%% Statistic
    bit_error_count = bit_error_count+sum(info~=info_hat);
    if mod(iter, 100)==0
        fprintf(repmat('\b', 1, 1+msg_length));
        msg = ['Compression Rate = ', num2str(rate), ', Run ', num2str(iter), ' times, BER = ', num2str(bit_error_count/N/iter, '%.2e')];
        msg_length = length(msg);
        fprintf([msg, '\n']);
    end
end


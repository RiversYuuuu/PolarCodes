function [send_signal] = polar_MLC_encoder(info, send_set, frozen_flag)
[level_num, N] = size(frozen_flag);
codeword = zeros(N, level_num);
for level = 1:level_num
    codeword(:, level) = polar_encoder_channel(info{1, level}, frozen_flag(level, :));
end
send_signal = send_set(1+bi2de(codeword)).';
end
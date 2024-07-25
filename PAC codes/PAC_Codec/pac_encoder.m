function [codeword] = pac_encoder(uncoded_bits, conv_vec)
N = length(uncoded_bits);
conv_length = length(conv_vec);
cache = zeros(1, conv_length);
conv_bits = zeros(1, N);
for index = 1:N
    cache = [uncoded_bits(index), cache(1:end-1)];
    conv_bits(index) = mod(sum(cache.*conv_vec), 2);
end
codeword = polar_encoder(conv_bits);
end
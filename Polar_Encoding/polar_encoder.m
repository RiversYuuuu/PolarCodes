function [codeword] = polar_encoder(uncoded_bits)
N = length(uncoded_bits);
n = log2(N);
index = (1:1:N);
odds = (1:2:N);
evens = (2:2:N);
uncoded_bits = reshape(uncoded_bits, N, 1);
temp = [uncoded_bits, uncoded_bits];
for i = 1:n
    odd_index = index(odds);
    even_index = index(evens);
    temp(even_index, 2) = temp(even_index, 1);
    temp(odd_index, 2) = mod(temp(even_index, 1)+temp(odd_index, 1), 2);
    index = reshape(index, 2, N/2);
    index = reshape(index',1, N);
    temp(:, 1) = temp(:, 2);
end
codeword = temp(:, 2);
end


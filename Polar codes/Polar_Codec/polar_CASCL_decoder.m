function [codeword_hat, info_hat] = polar_CASCL_decoder(LLR, list_size, CRCdetecter, type_flag, begin_layers, end_layers)

[codeword_set, uncoded_bits_set, ~] = polar_SCL_decoder4CRC_mex(LLR, list_size, type_flag, begin_layers, end_layers);

for list_index = 1:list_size
    info_tmp = uncoded_bits_set(type_flag==0, list_index);
    [info_tmp, error_flag] = CRCdetecter(info_tmp);
    if 1==list_index
        info_hat = info_tmp;
        codeword_hat = codeword_set(:, list_index);
    end
    if ~error_flag
        info_hat = info_tmp;
        codeword_hat = codeword_set(:, list_index);
        break;
    end
end

end
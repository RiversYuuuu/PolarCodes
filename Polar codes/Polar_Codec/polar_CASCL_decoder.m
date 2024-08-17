function [codeword_hat, info_hat] = polar_CASCL_decoder(LLR, list_size_range, CRCdetecter, type_flag, begin_layers, end_layers)
list_size = list_size_range(1);
while list_size<=list_size_range(end)
    % [codeword_set, uncoded_bits_set, uncoded_bits_soft] = polar_SCL_decoder4CRC(LLR, list_size, type_flag, begin_layers, end_layers);
    [codeword_set, uncoded_bits_set, uncoded_bits_soft] = polar_SCL_decoder4CRC_mex(LLR, list_size, type_flag, begin_layers, end_layers);
    for list_index = 1:list_size
        info_tmp = uncoded_bits_set(type_flag==0, list_index);
        [info_tmp, error_flag] = CRCdetecter(info_tmp);
        if 1==list_index || error_flag==0
            info_hat = info_tmp;
            codeword_hat = codeword_set(:, list_index);
            if ~error_flag
                break;
            end
        end
    end
    if ~error_flag
        break;
    else
        list_size = list_size*2;
    end
end
end

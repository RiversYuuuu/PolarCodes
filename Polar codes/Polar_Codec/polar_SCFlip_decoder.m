function [codeword_hat, info_hat] = polar_SCFlip_decoder(LLR, flip_time, frozen_flag, CRCdetector, begin_layers, end_layers)
[codeword_hat, uncoded_bits, uncoded_bits_soft] = polar_SC_decoder4flip(LLR, [], frozen_flag, begin_layers, end_layers);
info_hat = uncoded_bits(frozen_flag==0);
[info_hat, iserror] = CRCdetector(info_hat);
if iserror
    I_set = find(frozen_flag==0);
    I_soft = uncoded_bits_soft(frozen_flag==0);
    [~, sorted_indices] = sort(abs(I_soft));
    I_set = I_set(sorted_indices);
    for i = 1:flip_time(1)
        [codeword_tmp, uncoded_bits, ~] = polar_SC_decoder4flip(LLR, I_set(i), frozen_flag, begin_layers, end_layers);
        info_tmp = uncoded_bits(frozen_flag==0);
        [info_tmp, iserror] = CRCdetector(info_tmp);
        if ~iserror
            codeword_hat = codeword_tmp;
            info_hat = info_tmp;
            break;
        end
    end
%     if iserror
%         F_set = find(frozen_flag==1);
%         F_soft = uncoded_bits_soft(frozen_flag==1);
%         max_metric = sum(abs(uncoded_bits_soft(frozen_flag==0)));
%         max_info_hat = info_hat;
%         [~, sorted_indices] = sort(F_soft);
%         F_set = F_set(sorted_indices);
%         for i = 1:flip_time(2)
%             [codeword_hat, uncoded_bits, uncoded_bits_soft] = polar_SC_decoder4flip(LLR, F_set(i), frozen_flag, begin_layers, end_layers);
%             metric = sum(abs(uncoded_bits_soft(frozen_flag==0)));
%             info_hat = uncoded_bits(frozen_flag==0);
%             [info_hat, iserror] = CRCdetector(info_hat);
%             if ~iserror
%                 max_info_hat = info_hat;
%                 break;
%             end
%             if metric>max_metric
%                 max_metric = metric;
%                 max_info_hat = info_hat;
%             end
%         end
%         info_hat = max_info_hat;
%     end
end
end
function [mutual_info, mutual_info_var, capacity] = get_mutual_info(send_set, send_prob_set, SNR_range)
%%% Parameter
level_num = log2(length(send_set));
send_power = sum(send_prob_set.*abs(send_set).^2);
dim = 2-(sum(imag(send_set))==0);
mutual_info = zeros(level_num+1, length(SNR_range));

%%% Shannon Capacity
variance = 1./(dim*10.^(SNR_range/10));
capacity = dim*(1/2)*log2(1+1./variance);

%%% Quantization
quan_len = 0.05;
quan_max = 5*max(send_set);
quan_range = -quan_max:quan_len:quan_max;
recv_signal = zeros(1, length(quan_range)^dim);
recv_index = 1;
for recv_real = quan_range
    recv_signal(recv_index) = recv_real;
    for recv_imag = quan_range
        if dim==1
            recv_index = recv_index + 1;
            break
        end
        recv_signal(recv_index) = recv_real + recv_imag * 1i;
        recv_index = recv_index + 1;
    end
end

%%% P(X=0) Of Levels
levels_prob = cell(1, level_num);
for level = 1:level_num
    level_prob = zeros(2^level, 1);
    for index = 1:2^level
        step = 2^level;
        start_index = index;
        end_index = 2^level_num;
        level_prob(index) = sum(send_prob_set(start_index:step:end_index));
    end
    levels_prob{level} = level_prob;
end

%%% Calculate
for SNR_index = 1:length(SNR_range)
    SNR = SNR_range(SNR_index);
    AWGN_sigma = sqrt(send_power/(dim*10^(SNR/10)));
    %%% I(Y;bi|b1,b2,...,bi-1)=H(Y|b1,b2,...,bi-1)-H(Y|b1,b2,...,bi)
    for level = 1:level_num
        %%% H(Y|b1,b2,...,bi-1)
        H1 = 0;  
        for index = 1:2^(level-1)
            step = 2^(level-1);
            begin_index = index;
            end_index = 2^level_num;
            %%% P(Y|b1,b2,...,bi-1)
            signal_prob = reshape(send_prob_set(begin_index:step:end_index), [], 1);
            P1 = signal_prob.*get_AWGN_send_prob(send_set(begin_index:step:end_index), recv_signal, AWGN_sigma, dim);
            if size(P1, 1)>1
                P1 = sum(P1);
            end
            P1 = P1/sum(P1);
            %%% P(Y,b1,b2,...,bi-1)
            if level==1
                P2 = P1;
            else
                P2 = levels_prob{level-1}(index)*P1;
            end
            H1 = H1-sum(P2(P2>0).*log2(P1(P2>0)));
        end
        %%% H(Y|b1,b2,...,bi)
        H2 = 0;
        for index = 1:2^level
            step = 2^level;
            begin_index = index;
            end_index = 2^level_num;
            %%% P(Y|b1,b2,...,bi)
            signal_prob = reshape(send_prob_set(begin_index:step:end_index), [], 1);
            P3 = signal_prob.*get_AWGN_send_prob(send_set(begin_index:step:end_index), recv_signal, AWGN_sigma, dim);
            if size(P3, 1)>1
                P3 = sum(P3);
            end
            P3 = P3/sum(P3);
            %%% P(Y,b1,b2,...,bi)
            P4 = levels_prob{level}(index)*P3;
            H2 = H2-sum(P4(P4>0).*log2(P3(P4>0)));
        end
        mutual_info(level, SNR_index) = H1-H2;
    end
end
mutual_info(level_num+1,:) = sum(mutual_info(1:level_num,:));
mutual_info_var = var(mutual_info(1:level_num,:));
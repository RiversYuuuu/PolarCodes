function [send_prob] = get_AWGN_send_prob(send_set, recv_set, AWGN_sigma, dim)
    send_prob = zeros(length(send_set), length(recv_set));
    switch dim
        case 1
            for send_index=1:length(send_set)
                send_prob(send_index, :) = normpdf(send_set(send_index), recv_set, AWGN_sigma);
            end
        case 2
            for send_index=1:length(send_set)
                send_prob(send_index, :) = normpdf(real(send_set(send_index)), real(recv_set), AWGN_sigma);
                send_prob(send_index, :) = send_prob(send_index, :).*normpdf(imag(send_set(send_index)), imag(recv_set), AWGN_sigma);
            end
    end
end
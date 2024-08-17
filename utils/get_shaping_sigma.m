function [shaping_sigma, send_prob_set] = get_shaping_sigma(send_set)
dim = 2-(sum(abs(imag(send_set)))==0);
max_distance = max(max(real(send_set)), max(imag(real(send_set))));
shaping_sigma = ceil(max_distance/2.9/0.5)*0.5;
send_prob_set = reshape(get_AWGN_send_prob(send_set, 0, shaping_sigma, dim), 1, []);
send_prob_set = send_prob_set/sum(send_prob_set);
send_power = sum(send_prob_set.*abs(send_set).^2);
% fprintf('Desired power is %.2f, Actual power is %.2f\n', dim*shaping_sigma^2, send_power);
end
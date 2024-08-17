function [send_set, send_prob_set] = get_send_set(modulate_type, shaping_flag)
if isequal(modulate_type, '64QAM')
    send_set = qammod(0:63,64, get_SPmapping_order(64));
    send_prob_set = ones(1, 64)/64;
%     [mutual_info, mutual_info_var, capacity] = get_mutual_info(send_set, send_prob_set, 15);
end
if isequal(modulate_type, '64Spiral')
    m = 1:64;
    fs = 0;
    tm = sqrt(((4*pi*m).^2*fs/2)+sqrt((4*pi*m).^4*fs^2/4+(4*pi*m).^2));
    send_set = tm.*exp(tm*1i);
    if shaping_flag==0
%         order = get_maxVar_order(send_set, ones(1,64)/64, 10);
        order = [1, 2, 5, 4, 3, 6];
        index = de2bi(0:63);
        index = bi2de(index(:, order));
        send_set = send_set(1+index);
    end
    send_prob_set = ones(1, 64)/64;
end
end
function [P_up, P_down] = channel_polarization(P_array)
[~, channel_num]=size(P_array);
P_yox=zeros(2, channel_num);
P_yox(2, :) = P_array(1, :).*P_array(2, :);
P_yox(1, :) = P_array(1, :)-P_yox(2, :);
Pxy = kron(P_yox, P_yox);
I_up = [Pxy(1, :), Pxy(3, :); Pxy(4, :), Pxy(2, :)];
I_down = [Pxy(1, :)+Pxy(4, :); Pxy(2, :)+Pxy(3, :)];
P_up(1, :) = sum(I_up, 1);
P_up(2, :) = min(I_up)./sum(I_up, 1);
P_down(1, :) = sum(I_down, 1);
P_down(2, :) = min(I_down)./sum(I_down, 1);
end


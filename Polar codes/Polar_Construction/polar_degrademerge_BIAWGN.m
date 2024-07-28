function [error_bound, sorted_indices] = polar_degrademerge_BIAWGN(N, noise_sigma)
%%% Channel quantization
quan_level = 64;
quan_points = linspace(0, 1+2*noise_sigma, quan_level);
p0 = normpdf(quan_points, 1, noise_sigma);
p1 = normpdf(quan_points, -1, noise_sigma);
P_array = [(p0+p1)/sum(p0+p1); min([p0; p1])./(p0+p1)];
[error_bound, sorted_indices] = polar_degrademerge_BSC(N, P_array);
end
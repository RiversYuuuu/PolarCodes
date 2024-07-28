function [begin_layers, end_layers] = polar_decode_prepare(N)
n = log2(N);
begin_layers = ones(N, 1);
end_layers = zeros(N, 1);
begin_layers(1) = n;
for index = 2:N
    end_layer = 1;
    index_cp = index;
    while 1 == mod(index_cp,2)
        end_layer = end_layer + 1;
        index_cp = (index_cp+1)/2;
    end
    begin_layers(index,1) = end_layer;
    
    begin_layer = 0;
    index_cp = index;
    while 0 == mod(index_cp, 2)
        begin_layer = begin_layer + 1;
        index_cp = index_cp/2;
    end
    end_layers(index, 1) = begin_layer;
end
end


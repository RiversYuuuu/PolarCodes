function [entropy] = get_entropy(p)
    entropy = (1-p).*log2(1-p)+p.*(p);
end
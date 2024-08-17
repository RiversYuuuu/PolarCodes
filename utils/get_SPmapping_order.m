function order = get_SPmapping_order(M)
matrix = repmat("", sqrt(M), sqrt(M));
base = repmat(["0", "1"], 1, sqrt(M)/2);
base = [base, base(end:-1:1)];
for level = 1:log2(M)
    if 1==mod(level,2)
        matrix = reshape(matrix, 1, []);
        repete_num = 2^(round(level/2)-1);
        base_index = 0;
        for index = 1:repete_num:M
            matrix(index:index+repete_num-1) = base(1+base_index)+matrix(index:index+repete_num-1);
            base_index = mod(base_index+1, length(base));
        end
    else
        matrix = reshape(matrix, sqrt(M), sqrt(M));
        repete_num = 2^(level/2-1);
        flag = false;
        for row_index = 1:repete_num:sqrt(M)
            if flag==false
                matrix(row_index:row_index+repete_num-1,:) = "0"+matrix(row_index:row_index+repete_num-1,:);
            else
                matrix(row_index:row_index+repete_num-1,:) = "1"+matrix(row_index:row_index+repete_num-1,:);
            end
            flag = ~flag;
        end
    end
end
matrix = reshape(matrix, 1, []);
order = bin2dec(matrix);
end
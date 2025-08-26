function val = get_H(Y, test_type, c3) 
    if strcmp(test_type,'base')
        val = c3^2*(Y(:,1).^2 + Y(:,2).^2) - Y(:,3).^2;
    end
end


function val = get_H_der1(Y, dim_u, test_type, c3) %dim_u=1,2,3;
    if strcmp(test_type,'base')
        if dim_u == 1
            val = 2*c3^2*Y(:,1);
        elseif dim_u == 2
            val = 2*c3^2*Y(:,2);
        else
            val = -2*Y(:,3);
        end

    end
end


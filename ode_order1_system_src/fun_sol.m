function val = fun_sol(X, x_R, test_type, M, co_R, theta, deg_x, c3)  %only for x>0
    if strcmp(test_type,'base')
        len = length(X);
        val = zeros(len, 3);
        val(:,1) = X.^deg_x.*sin(theta*X);
        val(:,2) = X.^deg_x.*cos(theta*X);
        val(:,3) = c3*X.^deg_x;
        if len == M
            val(:,1) = val(:,1).*co_R;
            val(:,2) = val(:,2).*co_R;
            val(:,3) = val(:,3).*co_R;
        elseif len == 1 %find init with method = simple
            pos_temp = find(x_R ==X(1));
            val(:,1) = val(:,1).*co_R(pos_temp);
            val(:,2) = val(:,2).*co_R(pos_temp);
            val(:,3) = val(:,3).*co_R(pos_temp);
        end
    end
end


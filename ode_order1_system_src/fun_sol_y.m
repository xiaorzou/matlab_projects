function val = fun_sol_y(X, test_type,theta, deg_x, c3)  %only for x>0
    if strcmp(test_type,'base')
        len = length(X);
        val = zeros(len, 3);
        val(:,1) = X.^deg_x.*sin(theta*X);
        val(:,2) = X.^deg_x.*cos(theta*X);
        val(:,3) = c3*X.^deg_x;
    end
end


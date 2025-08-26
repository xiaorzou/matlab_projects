function val = fun_sol_rk(X, test_type, x_2N_R, co_2N_R, theta, deg_x, c3)  %only for x>0
    pos_temp = find(x_2N_R ==X(1));
    if strcmp(test_type,'base')
        len = length(X);
        val = zeros(len, 3);
        val(:,1) = X.^deg_x.*sin(theta*X)*co_2N_R(pos_temp);
        val(:,2) = X.^deg_x.*cos(theta*X)*co_2N_R(pos_temp);
        val(:,3) = c3*X.^deg_x*co_2N_R(pos_temp);
    end
end


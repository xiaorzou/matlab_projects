    function val = get_r_value_rk(X,x_2N_R, co_2N_R, co_2N_der_R, test_type, theta,deg_x, ps, qs, c3)
        len = length(X);
        pos_temp = find(x_2N_R ==X(1));
        if strcmp(test_type,'base')
            z_value = fun_der_sol_rk(X, test_type, x_2N_R, co_2N_R, co_2N_der_R, theta, deg_x, c3);
            u_value = fun_sol_rk(X, test_type, x_2N_R, co_2N_R, theta, deg_x, c3);
            val = zeros(len, 3);
            val(:,1) = z_value(:,1) - ps(1)* u_value(:,2).^2*co_2N_R(pos_temp) - qs(1)*u_value(:,1)*co_2N_R(pos_temp);
            val(:,2) = z_value(:,2) - ps(2)* u_value(:,3).^2*co_2N_R(pos_temp) - qs(2)*u_value(:,2)*co_2N_R(pos_temp);
            val(:,3) = z_value(:,3) - ps(3)* u_value(:,1).^2*co_2N_R(pos_temp) - qs(3)*u_value(:,3)*co_2N_R(pos_temp);
        end
    end


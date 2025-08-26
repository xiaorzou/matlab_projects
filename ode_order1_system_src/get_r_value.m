    function val = get_r_value(X,x_R, test_type, M, co_R, co_der_R,theta, deg_x, ps, qs, c3)
        len = length(X);
        if strcmp(test_type,'base')
            %z_value = fun_der_sol(X, M, co_R, co_der_R, co_2N_der_R, theta, deg_x, c3);
            z_value = fun_der_sol(X, x_R, M, co_R, co_der_R, theta, deg_x, c3);
            u_value = fun_sol(X,x_R, test_type, M, co_R, theta, deg_x, c3);
            val = zeros(len, 3);
            if len == M
                val(:,1) = z_value(:,1) - ps(1)*co_R.* u_value(:,2).^2 - qs(1)*co_R.*u_value(:,1);
                val(:,2) = z_value(:,2) - ps(2)*co_R.* u_value(:,3).^2 - qs(2)*co_R.*u_value(:,2);
                val(:,3) = z_value(:,3) - ps(3)*co_R.* u_value(:,1).^2 - qs(3)*co_R.*u_value(:,3);
            else
                val(:,1) = z_value(:,1) - ps(1)* u_value(:,2).^2 - qs(1)*u_value(:,1);
                val(:,2) = z_value(:,2) - ps(2)* u_value(:,3).^2 - qs(2)*u_value(:,2);
                val(:,3) = z_value(:,3) - ps(3)* u_value(:,1).^2 - qs(3)*u_value(:,3);
            end
        end
    end


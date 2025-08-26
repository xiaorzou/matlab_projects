    function val = fun_der_sol(X, x_R, M, co_R, co_der_R, theta, deg_x, c3)  %only for x>=0, direvative of u(x)*co_R
        len = length(X);
        val = zeros(len, 3);
        if deg_x>1
            val(:,1) = theta*X.^deg_x.*cos(theta*X) + (deg_x)*X.^(deg_x-1).*sin(theta*X) ;
            val(:,2) = -theta*X.^deg_x.*sin(theta*X) + (deg_x)*X.^(deg_x-1).*cos(theta*X);
            val(:,3) = c3*deg_x*X.^(deg_x-1);
        elseif deg_x ==0
            val(:,1) = theta*cos(theta*X);
            val(:,2) = -theta*sin(theta*X);
            val(:,3) = zeros(len,1);
        elseif deg_x ==1
            val(:,1) = theta*X.*cos(theta*X) + sin(theta*X);
            val(:,2) = -theta*X.*sin(theta*X) + cos(theta*X);
            val(:,3) = c3*ones(len,1);
        end
        
        if len==M
            val(:,1) = X.^deg_x.*sin(theta*X).*co_der_R + val(:,1).*co_R;
            val(:,2) = X.^deg_x.*cos(theta*X).*co_der_R + val(:,2).*co_R;
            val(:,3) = c3*X.^deg_x.*co_der_R + val(:,3).*co_R;
        elseif len == 1 %find init with method = simple
            pos_temp = find(x_R ==X(1));
            val(:,1) =  X.^deg_x.*sin(theta*X).*co_der_R(pos_temp) + val(:,1).*co_R(pos_temp);
            val(:,2) =  X.^deg_x.*cos(theta*X).*co_der_R(pos_temp) + val(:,2).*co_R(pos_temp);
            val(:,3) = c3*X.^deg_x.*co_der_R(pos_temp) + val(:,3).*co_R(pos_temp);
        end
    end


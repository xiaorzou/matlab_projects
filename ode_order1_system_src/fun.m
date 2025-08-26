    function val = fun(XU,x_R, test_type, M, co_R, co_der_R,theta, deg_x, ps, qs, c3)  % X>0 is must!  shape N , 4 X(:,1) x value,  X(:, 2) u1 value, X(:,3): u2 value, X(:,4): u3 value),  val shape: (V, 3) 
        [rows, temp] = size(XU);
        val = zeros(rows, 3);
        if strcmp(test_type,'base')
            
            r_value = get_r_value(XU(:,1),x_R, test_type, M, co_R, co_der_R,theta, deg_x, ps, qs, c3);
            if rows == M
                val(:,1) = r_value(:,1) + ps(1)* XU(:,3).^2.*co_R+qs(1)*XU(:,2).*co_R;
                val(:,2) = r_value(:,2) + ps(2)* XU(:,4).^2.*co_R+qs(2)*XU(:,3).*co_R;
                val(:,3) = r_value(:,3) + ps(3)* XU(:,2).^2.*co_R+qs(3)*XU(:,4).*co_R;            
            elseif rows == 1
                pos_temp = find(x_R ==XU(1));
                val(:,1) = r_value(:,1) + ps(1)* XU(:,3).^2*co_R(pos_temp)+qs(1)*XU(:,2)*co_R(pos_temp);
                val(:,2) = r_value(:,2) + ps(2)* XU(:,4).^2*co_R(pos_temp)+qs(2)*XU(:,3)*co_R(pos_temp);
                val(:,3) = r_value(:,3) + ps(3)* XU(:,2).^2*co_R(pos_temp)+qs(3)*XU(:,4)*co_R(pos_temp);        
            end
        end
    end


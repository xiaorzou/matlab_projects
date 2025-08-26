    function val = fun4rk(XU,x_2N_R, co_2N_R, co_2N_der_R, test_type, theta,deg_x, ps, qs, c3)  % X>0 is must!  shape N , 4 X(:,1) x value,  X(:, 2) u1 value, X(:,3): u2 value, X(:,4): u3 value),  val shape: (V, 3) 
        [rows, temp] = size(XU);
        val = zeros(rows, 3);
        pos_temp = find(x_2N_R ==XU(1));
        if strcmp(test_type,'base')
            r_value = get_r_value_rk(XU(:,1),x_2N_R, co_2N_R, co_2N_der_R, test_type, theta,deg_x, ps, qs, c3);
            val(:,1) = r_value(:,1) + ps(1)* XU(:,3).^2*co_2N_R(pos_temp)+qs(1)*XU(:,2)*co_2N_R(pos_temp);
            val(:,2) = r_value(:,2) + ps(2)* XU(:,4).^2*co_2N_R(pos_temp)+qs(2)*XU(:,3)*co_2N_R(pos_temp);
            val(:,3) = r_value(:,3) + ps(3)* XU(:,2).^2*co_2N_R(pos_temp)+qs(3)*XU(:,4)*co_2N_R(pos_temp);                           
        end
    end


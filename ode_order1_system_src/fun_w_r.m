    function val = fun_w_r(XU, r_value, test_type, ps, qs, M, co_R)  % X>0 is must!  shape N , 4 X(:,1) x value,  X(:, 2) u1 value, X(:,3): u2 value, X(:,4): u3 value),  val shape: (V, 3) 
        [rows, temp] = size(XU);
        val = zeros(rows, 3);
        if strcmp(test_type,'base')
            %r_value = get_r_value(XU(:,1));
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
        elseif strcmp(test_type,'rigid_body')
            if rows == M
                val(:,1) = ((1/I_rb(3) - 1/I_rb(2))* XU(:,3).*XU(:,4)).*co_R;
                val(:,2) = ((1/I_rb(1) - 1/I_rb(3))* XU(:,2).*XU(:,4)).*co_R;
                val(:,3) = ((1/I_rb(2) - 1/I_rb(1))* XU(:,2).*XU(:,3)).*co_R; 
            elseif rows == 1
                pos_this = find(x_R ==XU(1));
                val(:,1) = ((1/I_rb(3) - 1/I_rb(2))* XU(:,3).*XU(:,4))*co_R(pos_this);
                val(:,2) = ((1/I_rb(1) - 1/I_rb(3))* XU(:,2).*XU(:,4))*co_R(pos_this);
                val(:,3) = ((1/I_rb(2) - 1/I_rb(1))* XU(:,2).*XU(:,3))*co_R(pos_this);                    
            end
        end
    end

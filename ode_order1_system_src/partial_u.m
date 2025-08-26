function val = partial_u(XU, dim_F, dim_u, test_type, co_R, ps, qs)  % u_1' = R + (p_1*co_R) u_2^2 + (q_1*co_R) u_1
    len = length(XU(:,1));
    u1= XU(:,2);
    u2= XU(:,3);
    u3 = XU(:,4);
    if strcmp(test_type,'base')
        if dim_F == 1
            if dim_u == 1
                val = qs(1)*co_R.*ones(len,1);
            elseif dim_u == 2
                val = 2*ps(1)*co_R.*u2;
            else
                val = zeros(len,1);
            end
        elseif dim_F == 2
            if dim_u == 1
                val = zeros(len,1);
            elseif dim_u == 2
                val = qs(2)*co_R.*ones(len,1);
            else
                val = 2*ps(2)*co_R.*u3;
            end    
        else
            if dim_u == 1
                val = 2*ps(3)*co_R.*u1;
            elseif dim_u == 2
                val = zeros(len,1);
            else
                val = qs(3)*co_R.*ones(len,1);
            end  
        end
    end
    %val = val.*co_R;
end


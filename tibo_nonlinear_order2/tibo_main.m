
%%
%{
this is the driver to generate the numerical results shown in the following
paper

"Trigonometric Interpolation Based Optimization for Second Order Non-Linear
ODE  with Mixed Boundary Conditions"

Remark:
testing function and parameters
ODE
y'' = r(x) + p_uu y'^2 + p_uv*y'*y + p_vv y^2 + p_u y' + p_v y
y(x) = x^deg*cos(theta x)
deg = 0 and 1

there are five types of tests 
test_type = 1; 
    cover all four bondary conditions: bdy_types={'N', 'D', 'mix_1', 'mix_2'};
    one can set theta=pi/3 (para= -0.5+1/3),  or theta=3*pi/3 (para=1)

test_type = 2; 
    apply to bounary condition where y(s) is know,  and we have some
    knowedge about the range of y'(s).  
    
test_type = 3; 
    testing extra condition on v_u, and identify the solution with upper bound of v as 0.5
test_type = 4; %test extra condition on u
    testing extra condition on u_u, and
    identify the solution with upper bound of u as 0, i.e. deceasing over
    [s,e]
test_type = 5; %test extra condition on z
    testing extra condition on z_u, and
    identify the solution with upper bound of z as 0.5, i.e. concave downward over
    [s,e]
%}

function tibo_main()
    global cuf_off_para 
    cut_off_u = false;
    cut_off_v = false; 
    plot_only = true;   %default setting,  run general testing results for performance.
    plot_scen = 11;       %add flexibility to plot a single certain scenario 
    save_coef = false;
    save_v = false;
    save_u = false;
    save_z = false;
    condi_on_grid = true; %extra conditions are applied only on grid points over [s,e]. this is the default setting.  in practice, we only have knowledge of solution over [s,e]

    figtype_fig = '.fig';
    task = 'default';
     
    %test_type = 1; %default, for all tests except with extra conditions   
    %test_type = 2; %test extra condition on D_l/D_u
    %test_type = 3; %test extra condition on v
    %test_type = 4; %test extra condition on u
    test_type = 5; %test extra condition on z
    
    
    q = 7; %default setting in paper is 7
    deg = 1; 
    
    para = -0.5+1/3; 
    %para = 1; % two 
    
    theta = (0.5+para)*pi; 
    
    % initital setting of conditions,  don't change them!!!
    condi_Diri = 'no';  
    condi_v_low = false; 
    condi_v_up = false; 
    condi_u_low = false; 
    condi_u_up = false; 
    condi_z_low = false; 
    condi_z_up = false;   
    threshold_by_level = true;
    Diri_low_percent = 0;  % only applicable when condi_Diri = 'lowbound' with threshold_by_level = false;
    Diri_low_level = 0; % only applicable when condi_Diri = 'lowbound' with threshold_by_level = true;
    Diri_up_percent = 0;  % only applicable when condi_Diri = 'upbound' with threshold_by_level = false;
    Diri_up_level = 0;   % only applicable when condi_Diri = 'upbound' with threshold_by_level = true;
    condi_v_low_val = 0;
    condi_v_up_val = 0;
    condi_u_low_val = 0; 
    condi_u_up_val = 0;  
    condi_z_low_val = 0; 
    condi_z_up_val = 0;  
    low_val = 0;
    up_val = 0;
        
    
    if test_type == 1
        disp('default test, select proper parameters q, deg, para for impact analysis')
        extra_condi_type = 'none'; % extra condition default setting: extra_condi_type = 'none'
        %bdy_types = {'N', 'D', 'mix_1', 'mix_2'};% 'N' for 'Neumann' , D for 'Dirichlet', 'mix_1', 'mix_2';
        bdy_types = {'mix_2'};% 'N' for 'Neumann' , D for 'Dirichlet', 'mix_1', 'mix_2';
    elseif test_type >1  % test on other conditions        
        theta = pi/3;
        plot_scen = 7;
        if  test_type == 2
            bdy_types = {'D'};
            extra_condi_type = 'D_l';   
            if threshold_by_level
                Diri_low_level = -0.5;
            else
                Diri_low_percent= 0.10;
            end
        elseif test_type == 3
            bdy_types = {'mix_1'};
            extra_condi_type = 'v_u'; 
            condi_v_up_val = 0.5;
        elseif test_type == 4    
            bdy_types = {'mix_1'};
            extra_condi_type = 'u_u'; 
            condi_u_up_val = 0;
        elseif test_type == 5 
            bdy_types = {'mix_1'};
            extra_condi_type = 'z_u';
            condi_z_up_val = 0;
        end
    end
    
    if theta == 3*pi/2
        plot_location = 'northwest';
    elseif theta == pi/3
        plot_location = 'northeast';
    else
        plot_location = 'northwest';
    end                
    if strcmp(extra_condi_type,'D_u')
        condi_Diri = 'upbound';  
        up_val = Diri_up_percent;
    elseif strcmp(extra_condi_type,'D_l')
        condi_Diri = 'lowbound';  
        low_val = Diri_low_percent;
    elseif strcmp(extra_condi_type,'v_u')
        condi_v_up = true;    
        up_val = condi_v_up_val;
    elseif strcmp(extra_condi_type,'v_l')
        low_val = condi_v_low_val;
        condi_v_low = true;             
    elseif strcmp(extra_condi_type,'v_ul')
        condi_v_up = true;  
        condi_v_low = true;  
        low_val = condi_v_low_val;
        up_val = condi_v_up_val;
    elseif strcmp(extra_condi_type,'u_l')    
        low_val = condi_u_low_val;
        condi_u_low = true;     
    elseif strcmp(extra_condi_type,'u_u')        
        condi_u_up = true;   
        up_val = condi_u_up_val;
    elseif strcmp(extra_condi_type,'z_l')    
        low_val = condi_z_low_val;
        condi_z_low = true;         
    elseif strcmp(extra_condi_type,'z_u')        
        condi_z_up = true;     
        up_val = condi_z_up_val;
    end   

    fourier_normalizer_const()
    s = 1;
    e = 3*s;    
    p = q-1; %p should be smaller than q,  
    q_plot = 10; %10 default used in paper
    
    M = 2^q;
    N = 2*M;
    n = 2^p; % we should make n as large as possible,  as such, p=q-1 should be always true
    lambda = (e-s)/n;
    m = (M-n)/2;
    grid_pos = m+1:m+n+1;
    delta = lambda*m;
    o = s - delta;
    b = e-s + 2*delta;
    x_N = -b + o + lambda*(0:N-1);
    x_R = x_N(M+1:end);
    co_R = fourier_normalizer_cut_off(x_R, s-delta,s,e,e+delta,cuf_off_para); 
    p_para = {0.1, 0.1, 1, 0.1, 1}; % default   
    [p_uu, p_uv, p_vv, p_u, p_v] = p_para{:};
    
    %for cut-off of u and v
    b_u = 1000;
    b_v = 1000;
    del_v = 2;
    del_u = 2;
    l_v = -b_v;
    r_v = b_v;
    l_u = -b_u;
    r_u = b_u;
    s_v = l_v + del_v;
    e_v = r_v - del_v;
    s_u = l_u + del_u;
    e_u = r_u - del_u;
 
    
    %% testing specific functions functions 

   function val = get_sol_tibo(X)
        if strcmp(task, 'default')
            x = abs(X);
            val = x.^deg.*cos(theta*x);
        else
            disp(['not implemented with task ', task])
            return;
        end
    end

    function val = get_der1_sol_tibo(X)  
        if strcmp(task, 'default')
            x = abs(X);
            if deg==1
                val = cos(theta*x)-theta*x.*sin(theta*x);
            elseif deg==0
                val = -theta*sin(theta*x);
            else
                val = deg*x.^(deg-1).*cos(theta*x)-theta*x.^deg.*sin(theta*x);
            end
        else
            disp(['not implemented with task ', task])
            return;
        end
    end
    
    function val = get_f_tibo(X)
        if strcmp(task, 'default')
            x = X(1,:);
            v = X(2,:);
            u = X(3,:);
            if deg==1
                der_2 = -2*theta.*sin(theta*x)-theta^2*x.*cos(theta*x);
            elseif deg==2
                der_2 = deg*(deg-1)*cos(theta*x)-2*deg*theta*x.^(deg-1).*sin(theta*x)-theta^2*x.^deg.*cos(theta*x);
            elseif deg==0
                der_2 = -theta^2*cos(theta*x);
            else
                der_2 = deg*(deg-1)*x.^(deg-2).*cos(theta*x)-2*deg*theta*x.^(deg-1).*sin(theta*x)-theta^2*x.^deg.*cos(theta*x);
            end
            v_real = get_sol_tibo(x);
            u_real = get_der1_sol_tibo(x);
            if cut_off_v
                c_v = fourier_normalizer_cut_off(v, l_v,s_v,e_v,r_v,cuf_off_para);
            else
                c_v = ones(1,length(v));
            end
            if cut_off_u
                c_u = fourier_normalizer_cut_off(u, l_u,s_u,e_u,r_u,cuf_off_para);
            else
                c_u = ones(1,length(u));
            end
            val = der_2 - p_uu*u_real.^2 -p_uv*v_real.*u_real  - p_vv*v_real.^2 -p_v*v_real-p_u*u_real + p_uu*(c_u.*u).^2 + p_uv*c_v.*v.*c_u.*u + p_vv*(c_v.*v).^2 +  p_v*c_v.*v + p_u*(c_u.*u);
        else
            disp(['not implemented with task ', task])
            return;
        end
    end

    function val = get_r_tibo(x)
        if strcmp(task, 'default')
            if deg==1
                der_2 = -2*theta.*sin(theta*x)-theta^2*x.*cos(theta*x);
            elseif deg==2
                der_2 = deg*(deg-1)*cos(theta*x)-2*deg*theta*x.^(deg-1).*sin(theta*x)-theta^2*x.^deg.*cos(theta*x);
            elseif deg == 0
                der_2 = -theta^2*cos(theta*x);
            else
                der_2 = deg*(deg-1)*x.^(deg-2).*cos(theta*x)-2*deg*theta*x.^(deg-1).*sin(theta*x)-theta^2*x.^deg.*cos(theta*x);
            end
            v_real = get_sol_tibo(x);
            u_real = get_der1_sol_tibo(x);
            val = der_2 - p_uu*u_real.^2 -p_uv*v_real.*u_real  - p_vv*v_real.^2 -p_v*v_real-p_u*u_real;
        else
            disp(['not implemented with task ', task])
            return;
        end
    end
   
    function val = get_f_w_r_tibo(X, r)
        if strcmp(task, 'default')
            x = X(1,:);
            v = X(2,:);
            u = X(3,:);
            if cut_off_v
                c_v = fourier_normalizer_cut_off(v, l_v,s_v,e_v,r_v,cuf_off_para);
            else
                c_v = ones(1,length(v));
            end
            if cut_off_u
                c_u = fourier_normalizer_cut_off(u, l_u,s_u,e_u,r_u,cuf_off_para);
            else
                c_u = ones(1,length(u));
            end
            val = r + p_uu*(c_u.*u).^2 + p_uv*c_v.*v.*c_u.*u + p_vv*(c_v.*v).^2 +  p_v*c_v.*v + p_u*(c_u.*u);
        else
            disp(['not implemented with task ', task])
            return;
        end
    end

    function val = get_f_tibo_grid(X)   
        val = get_f_tibo(X);
        if length(val)==length(co_R)
            val = val.*co_R;
        else
            x1 = X(1,:);
            t=length(val);
            for i_t=1:t
                pos_this = find(x_R>=x1(i_t));
                val(i_t) =val(i_t)*co_R(pos_this(1));
            end
        end
    end

    function val = get_f_w_r_tibo_grid(X,r)   
        val = get_f_w_r_tibo(X,r);
        if length(val)==length(co_R)
            val = val.*co_R;
        else
            x1 = X(1,:);
            t=length(val);
            for i_t=1:t
                pos_this = find(x_R>=x1(i_t));
                val(i_t) =val(i_t)*co_R(pos_this(1));
            end
        end
    end
    
    function val = get_partial_u_tibo(X) %X=(x,v,u)
        if strcmp(task, 'default')
            v = X(2,:);
            u = X(3,:);
            if cut_off_u
                c_u = fourier_normalizer_cut_off(u, l_u,s_u,e_u,r_u,cuf_off_para);
                c_d_u = fourier_normalizer_cut_off_der1(u, l_u,s_u,e_u,r_u,cuf_off_para);
            else
                c_u = ones(1, length(u));
                c_d_u = zeros(1, length(u));
            end
            if cut_off_v
                c_v = fourier_normalizer_cut_off(v, l_v,s_v,e_v,r_v,cuf_off_para);
            else
                c_v = ones(1, length(v));
            end
            val = 2*p_uu*u.*c_u +2*p_uu*u.^2.*c_d_u + p_uv*c_v.*v.*u.*c_d_u + p_uv*c_v.*v.*c_u + p_u*c_u +p_u*u.*c_d_u;
           
        else
            disp(['not implemented with task ', task])
            return;
        end
    end

    function val = get_partial_u_tibo_grid(X)
        val = get_partial_u_tibo(X);
        if length(val)==length(co_R)
            val = val.*co_R;
        else
            t=length(val);
            x1 = X(1,:);
            for i=1:t
                pos_this = find(x_R>=x1(i));
                val(i) =val(i)*co_R(pos_this(1));
            end
        end
    end

    function val = get_partial_v_tibo(X) %X=(x,v,u)
        if strcmp(task, 'default')
            v = X(2,:);
            u = X(3,:);
            if cut_off_v
                c_v = fourier_normalizer_cut_off(v, l_v,s_v,e_v,r_v,cuf_off_para);
                c_d_v = fourier_normalizer_cut_off_der1(v, l_v,s_v,e_v,r_v,cuf_off_para);
            else
                c_v = ones(1, length(v));
                c_d_v = zeros(1, length(v));
            end
            
            %new
            if cut_off_u
               c_u = fourier_normalizer_cut_off(u, l_u,s_u,e_u,r_u,cuf_off_para);
            else
               c_u = ones(1, length(u));
            end
            val = 2*p_vv*v.*c_v + p_vv*v.^2.*c_d_v + p_uv*u.*c_u.*c_v + p_uv*u.*c_u.*v.*c_d_v + p_v*(c_v + v.*c_d_v);
        else
            disp(['not implemented with task ', task])
            return;
        end
    end

    function val = get_partial_v_tibo_grid(X)
        val = get_partial_v_tibo(X);
        if length(val)==length(co_R)
            val = val.*co_R;
        else
            t=length(val);
            x1 = X(1,:);
            for i=1:t
                pos_this = find(x_R>=x1(i));
                val(i) =val(i)*co_R(pos_this(1));
            end
        end
    end

%% start test
   
    z_1_true = get_sol_tibo(s);
    z_2_true = get_der1_sol_tibo(s);
    z_3_true = get_sol_tibo(e);
    z_4_true = get_der1_sol_tibo(e);
   
    cheat_derivative = get_der1_sol_tibo(s);
    cheat_fun = get_sol_tibo(s);
    if threshold_by_level
        if strcmp(condi_Diri,'lowbound')
            low_val = Diri_low_level;
        else
            up_val = Diri_up_level;
        end
    else
        if strcmp(condi_Diri,'lowbound')
            low_val = z_2_true - abs(z_2_true)*Diri_low_percent; %results applied in the paper!!!
        elseif strcmp(condi_Diri,'upbound')
            up_val = z_2_true + abs(z_2_true)*Diri_up_percent; 
        end
    end
    for jj = 1:length(bdy_types)        
        bdy_type = bdy_types(jj);
        scen_base = 1;
        scen_second = 1;
        if theta == pi/3
            if strcmp(bdy_type,'N')
                scen_base = 1; 
            elseif strcmp(bdy_type,'D')
                scen_second = 7; 
            elseif strcmp(bdy_type,'mix_1')
                scen_second = 14; 
            elseif strcmp(bdy_type,'mix_2')
                scen_base = 2;  
            end
        elseif theta == 3*pi/2
            if strcmp(bdy_type,'N')
                scen_base = 1; 
            elseif strcmp(bdy_type,'D')
                scen_second = 2; 
            elseif strcmp(bdy_type,'mix_1')
                scen_base = 4; 
            elseif strcmp(bdy_type,'mix_2')
                scen_second = 2;  
            end
        end
        rng('default')
        randoms = rand(2,5);
        %randoms = rand(2,1);
        init_guesses_der = randoms(1,:)-0.5;
        init_guesses_fun = randoms(2,:)-0.5;
        init_guesses_der = [cheat_derivative + init_guesses_der, ...
            2*cheat_derivative + init_guesses_der, -2*cheat_derivative + init_guesses_der,...
            3*cheat_derivative + init_guesses_der, -3*cheat_derivative + init_guesses_der];
        init_guesses_fun = [cheat_fun + init_guesses_fun, 2*cheat_fun + init_guesses_fun,...
            -2*cheat_fun + init_guesses_fun, 5*cheat_fun + init_guesses_fun, -5*cheat_fun + init_guesses_fun];

        scen_size = length(init_guesses_der);  
        output = zeros(scen_size, 16);
        if save_coef
            output_coef = zeros(M+1, scen_size);
        end
        if save_v
            output_v = zeros(n+1, scen_size);
        end
        if save_u
            output_u = zeros(n+1, scen_size);
        end
        if save_z
            output_z = zeros(n+1, scen_size);
        end        
        if strcmp(bdy_type, 'D')
            A1 = [1,0,0,0];
            A2 = [0,0,1,0]; 
        elseif strcmp(bdy_type, 'N')
            A1 = [1,0,0,0];
            A2 = [0,1,0,0]; 
        elseif strcmp(bdy_type, 'mix_1')
            A1 = [1,0,0,0];
            A2 = [0,0,0,1];
        elseif strcmp(bdy_type, 'mix_2')
            A1 = [1,1,0,0];
            A2 = [0,0,1,1];           
        end
        alpha = A1(1)*z_1_true+A1(2)*z_2_true+A1(3)*z_3_true+A1(4)*z_4_true;
        beta = A2(1)*z_1_true+A2(2)*z_2_true+A2(3)*z_3_true+A2(4)*z_4_true;
        AA = [A1;A2];
        greeks = [alpha, beta];
        for ii = 1:scen_size          
            if plot_only                            
                if ii~= plot_scen 
                    continue
                end            
            else
                if test_type ~= 1 && ii~=scen_base && ii~=scen_second
                    continue
                end
            end                            

            init_guess = [init_guesses_fun(ii),init_guesses_der(ii)];
            if strcmp(bdy_type, 'N')
                init_guess(1) = cheat_fun;
                init_guess(2) = cheat_derivative;
                if ii~=1
                    continue  %one scen is enough for N
                end
            end
            if strcmp(bdy_type, 'D')
                init_guess(1) = cheat_fun;
            end
            try
                [init_opt, myobjopt_shooting, X_rk4] = shooting_by_rk4(@get_f_tibo, s, e,  greeks, AA, n, init_guess); 
                [~, init_z] = runge_kutta_ode_order_k(@get_f_tibo_grid, x_R,init_opt, m+1, []); 
            catch
                warning('failing in shooting_by_rk4 in scen %d.  using zeros as inti in opt', ii);
                X_rk4 = zeros(2,M/2+1);
                %init_opt = zeros(1,2);
                myobjopt_shooting = -1;
                init_z = zeros(1,M);
            end
            
            if isinf(max(abs(init_z))) || test_type==0
                init_z = zeros(1,M);
            end                            
            r_values = get_r_tibo(x_R);
            [X_inte, a_0, a_1, b,o,myobjopt_inte] = tibo_engine(@get_f_w_r_tibo_grid, ...
                @get_partial_u_tibo_grid, @get_partial_v_tibo_grid,...
                s, e, p,q, greeks, AA,init_z, r_values,...
                condi_on_grid, condi_Diri, condi_v_low, condi_v_up, condi_u_low,...
                condi_u_up, condi_z_low,  condi_z_up, up_val, low_val);
            if myobjopt_inte == -1
                continue
            end
            
            if save_v
                output_v(:,ii) = X_inte(2,grid_pos)';
            end
            if save_u
                output_u(:,ii) = X_inte(3,grid_pos)';
            end
            if save_z
                output_z(:,ii) = X_inte(4,grid_pos)';
            end            
            
            z_N = X_inte(4,:);
            z_N = [0,-fliplr(z_N(2:M)), z_N];  
            J_M = zeros(1,M);
            J_M(2:M) = 1./(1:M-1);
            J_N = zeros(1,N);
            J_N(1:M) = J_M;
            A = ones(1,N);
            A(2:2:end) = -1;
            b_M = 2*imag(ifft(z_N));
            b_M = b_M.*A;
            b_M = b_M(1:M);
            bopi = b/pi;
            coef_u = zeros(1,M);
            coef_u(1) = a_0;
            coef_u(2:end) = -bopi*J_N(2:M).*b_M(2:M);
            coef_v = zeros(1,M);
            coef_v(2:M) = -bopi^2*J_N(2:M).^2.*b_M(2:M);

            x_plot = X_inte(1,:);
            pos_temp = find(x_plot<=e & x_plot>=s);
            x_plot = x_plot(pos_temp);
            v_plot = X_inte(2,:);
            u_plot = X_inte(3,:);
            v_plot = v_plot(pos_temp);
            u_plot = u_plot(pos_temp);
            v_true = get_sol_tibo(x_plot);
            v_error_inte = v_plot - v_true;
            v_error_rk4 = X_rk4(1,:) - v_true;            
            alpha_inte = A1(1)*v_plot(1)+A1(2)*u_plot(1)+A1(3)*v_plot(end)+A1(4)*u_plot(end);
            beta_inte = A2(1)*v_plot(1)+A2(2)*u_plot(1)+A2(3)*v_plot(end)+A2(4)*u_plot(end);                      
            if save_coef
                if abs(myobjopt_inte)< 10^-6
                    output_coef(1:2, ii) = [a_1, a_0];
                    output_coef(3:end, ii) = -bopi^2*J_N(2:M).^2.*b_M(2:M);
                end
            end
            %q_plot = q;
            p_plot = q_plot-1; 
            M_plot = 2^q_plot;
            N_plot = 2*M_plot;
            n_plot = 2^p_plot; 
            lambda_plot = (e-s)/n_plot;
            x_N_plot = -b + o + lambda_plot*(0:N_plot-1);
            x_R_plot = x_N_plot(M_plot+1:end);
            u_M_plot = cos_approx_engine_coef2value_general(coef_u, x_R_plot, b, 'cos');
            v_M_plot = cos_approx_engine_coef2value_general(coef_v, x_R_plot, b, 'sin');
            v_M_plot = v_M_plot + a_1 + a_0.*x_R_plot;
            z_left_plot = cos_approx_engine_coef2value_general(b_M, x_R_plot, b, 'sin');

            X = zeros(3, M_plot);
            X(1,:) = x_R_plot;
            X(2,:) = v_M_plot;
            X(3,:) = u_M_plot;
            z_right_plot = get_f_tibo(X);
            diff = z_left_plot-z_right_plot;
            pos_temp = find(x_R_plot<=e & x_R_plot>=s);
            diff = diff(pos_temp);
            if ii == scen_base || ii==scen_second || ii == plot_scen
                v_true_plot = get_sol_tibo(x_R_plot);
                xlabel_this = '$x$';   
                filename_fig = strcat('output/tibo_', bdy_type, '_type_', num2str(test_type),'_', ...
                    extra_condi_type, '_', num2str(low_val),...
                    '_' , num2str(up_val),  '_',  num2str(condi_on_grid), '_', num2str(ii),...
                    '_', num2str(q), '_', num2str(deg),...
                    '_', num2str(round(theta)), figtype_fig);
                title = strcat('The $v$ and $y$ over $[0,b]$');
                legend_y = '$v$';
                legend_z = '$y$';
                ylabel_this = '';
                plot_latex_2(filename_fig{1},x_R_plot, v_M_plot,  v_true_plot, title,...
                    legend_y, legend_z,  xlabel_this, ylabel_this, plot_location);
            end
            output(ii,:) = [ii, q, theta, deg, p_uu, p_uv, p_vv, p_u, p_v, ...
                max(abs(v_error_inte)), max(abs(v_error_rk4)), max(abs(diff)),...
                myobjopt_inte, myobjopt_shooting, alpha_inte-alpha, beta_inte-beta];
        end
        if (~plot_only)
            vnames = {'scen_id', 'q', 'theta', 'deg', 'p_uu', 'p_uv', 'p_vv', 'p_u', 'p_v',...
                'v_err_inte', 'v_err_rk4', 'max_diff', ...
                'fval_opt', 'fval_shooting', 'err_alpha', 'err_beta'};
            tibo_o_file = strcat('output/tibo_', bdy_type, '_type_', num2str(test_type), '_c_', extra_condi_type, '_l_', num2str(low_val),...
                        '_u_' , num2str(up_val),  '_c_',  num2str(condi_on_grid), '_q_', num2str(q), '_d_', num2str(deg), '_t_', num2str(round(theta)),  '.xlsx');
            mylib_writearray(vnames, output, tibo_o_file{1});
            vnames = {};
            for iii = 1:scen_size
                vnames = [vnames, strcat('scen_',num2str(iii))];
            end
            if save_coef
                output_file = strcat('output/tibo_coef_', bdy_type, '_type_', num2str(test_type), '_c_', extra_condi_type, '_l_', num2str(low_val),...
                            '_u_' , num2str(up_val),  '_c_',  num2str(condi_on_grid), '_q_', num2str(q), '_d_', num2str(deg), '_t_', num2str(round(theta)), '.xlsx');
                mylib_writearray(vnames, output_coef, output_file{1});
            end
            if save_v
                output_file = strcat('output/tibo_v_', bdy_type, '_type_', num2str(test_type), '_c_', extra_condi_type, '_l_', num2str(low_val),...
                            '_u_' , num2str(up_val),  '_c_',  num2str(condi_on_grid), '_q_', num2str(q), '_d_', num2str(deg), '_t_', num2str(round(theta)), '.xlsx');
                mylib_writearray(vnames, output_v, output_file{1});            
            end
            if save_u
                output_file = strcat('output/tibo_u_', bdy_type, '_type_', num2str(test_type), '_c_', extra_condi_type, '_l_', num2str(low_val),...
                            '_u_' , num2str(up_val),  '_c_',  num2str(condi_on_grid), '_q_', num2str(q), '_d_', num2str(deg), '_t_', num2str(round(theta)), '.xlsx');
                mylib_writearray(vnames, output_u, output_file{1});            
            end
            if save_z
                output_file = strcat('output/tibo_z_', bdy_type, '_type_', num2str(test_type), '_c_', extra_condi_type, '_l_', num2str(low_val),...
                            '_u_' , num2str(up_val),  '_c_',  num2str(condi_on_grid), '_q_', num2str(q), '_d_', num2str(deg), '_t_', num2str(round(theta)), '.xlsx');
                mylib_writearray(vnames, output_z, output_file{1});            
            end                                               
        end
    end
end
    

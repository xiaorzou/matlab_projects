
%{

testing function and parameters

ODE
y'' = r(x) + 

y(x) = x^deg*cos(theta x)
deg = 0 and 1
theta=pi/2, 3*pi/2

tset type:

%}

function fourier_normalizer_app_ode_order2_en_test()
    %addpath('D:/matlab/fourier_normalizer/fourier_normalizer_engine')  
    %addpath('D:/matlab/cos_approx/cos_approx_engine')  
    %addpath('D:/matlab/mylib') 
    global cuf_off_para 
    %cut_off_u = true;
    %cut_off_v = true;
    cut_off_u = false;
    cut_off_v = false;
    flag_display = false;
    %flag_plot = true;
    init_by_neumann = false;
    plot_only = false;
    
    %apply_condition_diri = 'upbound';
    apply_condition_diri = 'no';
    %apply_condition_diri = 'lowbound'; %'no', 'upbound'
    apply_condition_on_lower_boundary_flag = false;
    apply_condition_on_lower_boundary_value = -0.1;
    apply_condition_on_upper_boundary_flag = false;
    apply_condition_on_upper_boundary_value = 10^6;

    
    figtype_fig = '.fig';
    %task = 'linear';
    task = 'non_linear';
    %init_method = 'runge_kutta';
    plot_location = 'northwest';
    constrainttolerance = 10^-15;
    maxiter = 1000;

    fourier_normalizer_const()
    s = 1;
    e = 3*s;
    q = 4; %default setting in paper
    %q = 10;
    p = q-1; %p should be smaller than q,  
    q_plot = 10; %default used in paper
    
    deg = 1; %donot change it.  deg=2 might not work
    para = 1; %or para =1
    theta = (0.5+para)*pi;
    
    plot_scen = 7; %default used in paper
    
    M = 2^q;
    N = 2*M;
    n = 2^p; % we should make n as large as possible,  as such, p=q-1 should be always true
    lambda = (e-s)/n;
    m = (M-n)/2;
    delta = lambda*m;
    o = s - delta;
    b = e-s + 2*delta;
    x_N = -b + o + lambda*(0:N-1);
    x_R = x_N(M+1:end);
     
    co_R = fourier_normalizer_cut_off(x_R, s-delta,s,e,e+delta,cuf_off_para);
    
    %test_case = 1; % v.^2 + 0.1*u.^2;   temperory!!! only used to test non-uniquesnss! also,
    %test_case = 2; % v + 0.1*u.^2; keep it! 
    %test_case = 3; % v.^2 + 0.1*u;  
    %test_case = 4; % 0.1*v.*u;   
    %test_case = 5; % p_vv*v^2 + p_v v + p_u u;   
    %para = 0 is must, otherwise,  inf will appear in initial values. 
    
    %para_scen = 'v2uv';
    para_scen = 'u2v2uvuv';
    %para_scen = 'u2uvuv_neg_v';
    %para_scen = 'u2uv';
    %para_scen = 'uvuv';
    %para_scen = 'vv_small';
    %para_scen = 'vv';
    %para_scen = 'uv';  %test linear ODE
    %para_scen = 'u2uvuv';
    %cond_types = {'Neumann','Dirichlet', 'mix_1', 'mix_2'};
    %cond_types = {'Neumann','Dirichlet'};
    %cond_types = {'Neumann'};
    %cond_types = {'Dirichlet'};
    cond_types = {'mix_1'};
    %cond_types = {'mix_2'};
    
    
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
    

    
    v_max = fun_sol(x_R); 
    v_max = max(abs(v_max));

    
    if strcmp(para_scen, 'v2uv')
        p_para = {0, 0, 1, 0.1, 1}; % default
    elseif strcmp(para_scen, 'u2v2uvuv')
        p_para = {0.1, 0.1, 1, 0.1, 1}; % default        
    elseif strcmp(para_scen, 'u2uv')
        p_para = {0.1, 0, 0, 0.1, 1}; % default        
    elseif strcmp(para_scen, 'uvuv')
        p_para = {0, 1, 0, 0.1, 1}; % default
    elseif strcmp(para_scen, 'vv_small')
        p_para = {0.0, 0.0, 6/(v_max*4*b^2/2), 0.0, 0.0}; % not unique,  if /2 is removed, unique
    elseif strcmp(para_scen, 'vv')
        p_para = {0.0, 0.0, 1, 0, 0}; % default    
    elseif strcmp(para_scen, 'uv')
        p_para = {0.0, 0.0, 0, 0.1, 1}; % default  
    elseif strcmp(para_scen, 'u2uvuv')
        p_para = {0.1, 0.1, 0, 0.1, 1}; % default  
    elseif strcmp(para_scen, 'u2uvuv_neg_v')    
        p_para = {0.1, 0.1, 0, 0.1, -1}; % default  
    end
    [p_uu, p_uv, p_vv, p_u, p_v] = p_para{:};
    function val = fun_sol(X)    %only for x<=0 old:   change to x>=0
        if strcmp(task, 'non_linear')
            x = abs(X);
            val = x.^deg.*cos(theta*x);
        else
            disp(['not implemented with task ', task])
            return;
        end
    end

    function val = fun_der_sol(X)  
        if strcmp(task, 'non_linear')
            x = abs(X);
            if deg==1
                val = cos(theta*x)-theta*x.*sin(theta*x);
            else
                val = deg*x.^(deg-1).*cos(theta*x)-theta*x.^deg.*sin(theta*x);
            end
        else
            disp(['not implemented with task ', task])
            return;
        end
    end

    function val = fun_der2_sol(X)  %only for x<=0
        if strcmp(task, 'non_linear')
            x = abs(X);
            if deg==1
                val = -2*theta*sin(theta*x)-theta^2*x.*cos(theta*x);
            else
                val = deg*(deg-1)*x.^(deg-2).*cos(theta*x)-2*deg*theta*x.^(deg-1).*sin(theta*x)-theta^2*x.^deg.*cos(theta*x);
            end
        else
            disp(['not implemented with task ', task])
            return;
        end
    end

    %F(x,y',y):F(x,u,v)
    %F(x,u,v) = deg*(deg-1)*x.^(deg-2).*cos(theta*x)-2*deg*theta*x.^(deg-1).*sin(theta*x)-theta^2*x.^deg.*cos(theta*x) - (x^d cos(theta x))^2 - (deg*x.^(deg-1).*cos(theta*x)-theta*x.^deg.*sin(theta*x))^2 + 0.1u^2 + v^2;
    
    function val = fun(X)
        if strcmp(task, 'non_linear')
            x = X(1,:);
            v = X(2,:);
            u = X(3,:);
            if deg==1
                der_2 = -2*theta.*sin(theta*x)-theta^2*x.*cos(theta*x);
            elseif deg==2
                der_2 = deg*(deg-1)*cos(theta*x)-2*deg*theta*x.^(deg-1).*sin(theta*x)-theta^2*x.^deg.*cos(theta*x);
            else
                der_2 = deg*(deg-1)*x.^(deg-2).*cos(theta*x)-2*deg*theta*x.^(deg-1).*sin(theta*x)-theta^2*x.^deg.*cos(theta*x);
            end
            v_real = fun_sol(x);
            u_real = fun_der_sol(x);
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
            %{
            if test_case == 1
                val = der_2 - p_vv*v_real.^2 - p_uu*u_real.^2  + p_vv*(c_v.*v).^2 + p_uu*(c_u.*u).^2; 
            elseif test_case == 2
                val = der_2 - p_v*v_real - p_uu*u_real.^2  + p_v*c_v.*v + p_uu*(c_u.*u).^2;
            elseif test_case == 3
                val = der_2 - p_vv*v_real.^2 - p_u*u_real  + p_vv*c_v.*v.^2 + p_u*(c_u.*u);
            elseif test_case == 4
                val = der_2 - p_uv*v_real.*u_real  + p_uv*c_v.*v.*c_u.*u;
            elseif test_case == 5
                val = der_2 - p_vv*v_real.^2 -p_v*v_real-p_u*u_real + p_vv*(c_v.*v).^2 +  p_v*c_v.*v + p_u*(c_u.*u);
            end
            %}
        else
            disp(['not implemented with task ', task])
            return;
        end
    end

    function val = fun_grid(X)   
        val = fun(X);
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

    function val = partial_u(X) %X=(x,v,u)
        if strcmp(task, 'non_linear')
            v = X(2,:);
            u = X(3,:);
            if cut_off_u
                c_u = fourier_normalizer_cut_off(u, l_u,s_u,e_u,r_u,cuf_off_para);
                c_d_u = fourier_normalizer_cut_off_der1(u, l_u,s_u,e_u,r_u,cuf_off_para);
            else
                c_u = ones(1, length(u));
                c_d_u = zeros(1, length(u));
            end
            
            %new 
            if cut_off_v
                c_v = fourier_normalizer_cut_off(v, l_v,s_v,e_v,r_v,cuf_off_para);
            else
                c_v = ones(1, length(v));
            end
            val = 2*p_uu*u.*c_u +2*p_uu*u.^2.*c_d_u + p_uv*c_v.*v.*u.*c_d_u + p_uv*c_v.*v.*c_u + p_u*c_u +p_u*u.*c_d_u;
            
            %val = 0.2*u;
            %{
            if test_case == 1 || test_case == 2
                val = 2*p_uu*u.*c_u +2*p_uu*u.^2.*c_d_u;
            elseif test_case == 3 || test_case == 5
                val = p_u*c_u +p_u*u.*c_d_u;
            elseif test_case == 4
                v = X(2,:);
                if cut_off_v
                    c_v = fourier_normalizer_cut_off(v, l_v,s_v,e_v,r_v,cuf_off_para);
                else
                    c_v = ones(1, length(v));
                end
                %val = 0.1*c_v.*v.*u.*c_d_u + 0.1*c_v.*v.*c_u;
                val = p_uv*c_v.*v.*u.*c_d_u + p_uv*c_v.*v.*c_u;
            else
                disp(['invalid test case ', test_case]);
                return
            end
            %}
        else
            disp(['not implemented with task ', task])
            return;
        end
    end

    function val = partial_v(X) %X=(x,v,u)
        if strcmp(task, 'non_linear')
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
            %{
            if test_case == 1 || test_case == 3
                %val = 2*v
                %val = 2*v.*c_v + v.^2.*c_d_v;
                val = 2*p_vv*v.*c_v + p_vv*v.^2.*c_d_v;
            elseif test_case == 2
                %val = 1;
                %val = c_v + v.*c_d_v;
                val = p_v*(c_v + v.*c_d_v);
            elseif test_case == 4
                u = X(3,:);
                if cut_off_u
                    c_u = fourier_normalizer_cut_off(u, l_u,s_u,e_u,r_u,cuf_off_para);
                else
                    c_u = ones(1, length(u));
                end
                %val = 0.1*u.*c_u.*c_v + 0.1*u.*c_u.*v.*c_d_v;
                val = p_uv*u.*c_u.*c_v + p_uv*u.*c_u.*v.*c_d_v;
            elseif test_case == 5
                %val = 0.1*u.*c_u.*c_v + 0.1*u.*c_u.*v.*c_d_v;
                val = 2*p_vv*v.*c_v + p_vv*v.^2.*c_d_v + p_v*(c_v + v.*c_d_v);
            end
            %}
        else
            disp(['not implemented with task ', task])
            return;
        end
    end
    % F = g(x) + p(x)*v^2 + q(x)*v + r(x)*u  %not needed
    % F = g(x) + q(x)*v + r(x)*u  % this is tested.  linear 
    deg_p = 0;
    deg_q = 0;
    deg_r = 0;
    coef_p = 0;
    coef_q = 1;
    coef_r = 1;
    
    function val = generate_g(x_R)
        z =  fun_der2_sol(x_R);
        v = fun_sol(x_R);
        u = fun_der_sol(x_R);
        %val = z-p_vv*v.^2-p_v*v-p_u*u;
        val = z-p_v*v-p_u*u; %only linear term!
    end
    %fun_p is not needed
    function val = fun_p(X)
        %val = p_vv*ones(1,length(X)).*co_R;
        val = p_vv*ones(1,length(X));
    end
    function val = fun_q(X)
        %val = p_v*ones(1,length(X)).*co_R;
        val = p_v*ones(1,length(X));
    end

    function val = fun_r(X)
        %val = p_u*ones(1,length(X)).*co_R;
        val = p_u*ones(1,length(X));
    end
    z_1_true = fun_sol(s);
    z_2_true = fun_der_sol(s);
    z_3_true = fun_sol(e);
    z_4_true = fun_der_sol(e);
   
    cheat_derivative = fun_der_sol(s);
    cheat_fun = fun_sol(s);
    if strcmp(apply_condition_diri,'lowbound')
        cond_b = z_2_true - abs(z_2_true)*0.10;
    elseif strcmp(apply_condition_diri,'upbound')
        cond_b = z_2_true + abs(z_2_true)*0.25; 
    else
        cond_b = 0;
    end
    for jj = 1:length(cond_types)
        cond_type = cond_types(jj);
        rng('default')
        randoms = rand(2,5);
        %randoms = rand(2,1);
        init_guesses_der = randoms(1,:)-0.5;
        init_guesses_fun = randoms(2,:)-0.5;
        init_guesses_der = [cheat_derivative + init_guesses_der, 2*cheat_derivative + init_guesses_der, -2*cheat_derivative + init_guesses_der, 3*cheat_derivative + init_guesses_der, -3*cheat_derivative + init_guesses_der];
        init_guesses_fun = [cheat_fun + init_guesses_fun, 2*cheat_fun + init_guesses_fun, -2*cheat_fun + init_guesses_fun, 5*cheat_fun + init_guesses_fun, -5*cheat_fun + init_guesses_fun];

        output = zeros(length(init_guesses_der), 12);
        output_max_error = zeros(length(init_guesses_der), 8);
        output_coef = zeros(M+1, length(init_guesses_der));
        if strcmp(cond_type, 'Dirichlet')
            % 'Dirichlet'
            A1 = [1,0,0,0];
            A2 = [0,0,1,0]; 
            %alpha = fun_sol(s);
            %beta = fun_sol(e);
        elseif strcmp(cond_type, 'Neumann')
            A1 = [1,0,0,0];
            A2 = [0,1,0,0]; 
            %alpha = fun_sol(s);
            %beta = fun_der_sol(s);
        elseif strcmp(cond_type, 'mix_1')
            %non-linear
            %A1 = [1,1,0,0];
            %A2 = [0,0,1,1];
            %linear
            A1 = [1,0,0,0];
            A2 = [0,0,0,1];
            
            
            %alpha = fun_sol(s);
            %beta = fun_der_sol(e);
        elseif strcmp(cond_type, 'mix_2')
            %nonlinear
            %A1 = [1,0,1,0]; 
            %A2 = [0,1,0,1];
            %linear
            A1 = [1,1,0,0];
            A2 = [0,0,1,1];
        elseif strcmp(cond_type, 'mix_3')
            A1 = [1,0,0,0]; 
            A2 = [0,0,0,1];
        end
        alpha = A1(1)*z_1_true+A1(2)*z_2_true+A1(3)*z_3_true+A1(4)*z_4_true;
        beta = A2(1)*z_1_true+A2(2)*z_2_true+A2(3)*z_3_true+A2(4)*z_4_true;
        AA = [A1;A2];
        greeks = [alpha, beta];
        for ii = 1:length(init_guesses_der)
            sprintf('scen: %d\n', ii)
            
            if ii~= plot_scen && plot_only
                continue
            end
            if ii~= 11
                %continue
            end
            init_guess = [init_guesses_fun(ii),init_guesses_der(ii)];
            if strcmp(cond_type, 'Neumann')
                init_guess(1) = cheat_fun;
                init_guess(2) = cheat_derivative;
                if ii~=1
                    continue
                end
            end
            if strcmp(cond_type, 'Dirichlet')
                init_guess(1) = cheat_fun;
            end
            
            %try
            
            %[init_opt, myobjopt_shooting, X_rk4] = fourier_normalizer_app_ode_get_init_value_by_rk4_order2(@fun, s, e,  greeks, AA, n, init_guess, delta, maxiter, constrainttolerance);  
            [init_opt, myobjopt_shooting, X_rk4] = fourier_normalizer_app_ode_get_init_value_by_rk4_order2(@fun, s, e,  greeks, AA, n, init_guess, maxiter, constrainttolerance); 
            [~, init_z] = runge_kutta_ode_order_k(@fun_grid, x_R,init_opt, m+1, []); 
                if init_by_neumann
                    AA_n  = [1,0,0,0;0,1,0,0];                 
                    [X_inte_n,~,~,~,~,~] = fourier_normalizer_app_ode_order2_en(@fun, @partial_u, @partial_v, s, e, p,q, init_opt, AA_n,...
                        flag_display,  constrainttolerance,maxiter, apply_condition_diri,cond_b,...
                        apply_condition_on_lower_boundary_flag, apply_condition_on_lower_boundary_value, apply_condition_on_upper_boundary_flag, apply_condition_on_upper_boundary_value, init_z);
                    init_z = X_inte_n(4,:);
                end
            
                [X_inte, a_0, a_1, b,o,myobjopt_inte] = fourier_normalizer_app_ode_order2_en(@fun, @partial_u, @partial_v, s, e, p,q, greeks, AA, ...
                    flag_display,  constrainttolerance,maxiter, apply_condition_diri,cond_b,...
                    apply_condition_on_lower_boundary_flag, apply_condition_on_lower_boundary_value, apply_condition_on_upper_boundary_flag, apply_condition_on_upper_boundary_value, init_z);
                
                x_M = X_inte(1,:);
                %x_R_plusb = zeros(1, M+1);
                %x_R_plusb(1:M) = x_M;
                %x_R_plusb(M+1) = b;
                %v_R_plusb_true = fun_sol(x_R_plusb);
                
                %P = fun_p(x_M);
                Q = fun_q(x_M);
                Q = Q.*co_R;
                R = fun_r(x_M);
                R = R.*co_R;
                %G = X_inte(4,:)-(X_inte(2,:).^2).*P - X_inte(2,:).*Q -  X_inte(3,:).*R;
                
                
                %G = X_inte(4,:)- X_inte(2,:).*Q -  X_inte(3,:).*R;
                %G = fun_der2_sol(x_M).*co_R - X_inte(2,:).*Q -  X_inte(3,:).*R;
                %G_check = generate_g(x_M);
                G = generate_g(x_M);
                G = G.*co_R;
                %pos_temp = find(x_R>=s & x_R<=e);
                %disp(max(abs(G(pos_temp)-G_check(pos_temp))));
                %if strcmp(cond_type, 'Dirichlet') && strcmp(para_scen, 'uv')
                max_error_linear = 0;
                if  strcmp(para_scen, 'uv') && abs(myobjopt_inte)<10^(-6)
                    %V_linear = investigate(X_inte, a_0, a_1, b, G, P, Q, R,  alpha, beta, m, n, cond_type,  A1,A2);
                    V_linear = investigate(X_inte, a_0, a_1, b, G, Q, R,  alpha, beta, m, n, cond_type,  A1,A2);
                    x_this = X_inte(1,:);    
                    x_R_plusb = zeros(1, M+1);
                    x_R_plusb(1:M) = x_this;
                    x_R_plusb(M+1) = b;
                    v_true_this = fun_sol(x_R_plusb);
                    v_Mplus1 = zeros(1,M+1);
                    v_Mplus1(1:M) = X_inte(2,:);
                    v_Mplus1(M+1) = a_1 + a_0*b;
                    pos_this_ = find(x_R_plusb>=s & x_R_plusb<=e);
                    max_error_linear = max(abs(v_true_this(pos_this_)'-V_linear(pos_this_)));
                    max_error_non_linear = max(abs(v_true_this(pos_this_)'-v_Mplus1(pos_this_)'));
                    res_inv = max(abs(V_linear(pos_this_)'-v_Mplus1(pos_this_)));
                    fprintf('max_error_linear %e, max_error_non_linear %e  max_compare %e.\n',max_error_linear,max_error_non_linear, res_inv);
                    %{
                    uv test,  scen 1
                    Neauman: max_error_linear 1.0976e-09, max_error_non_linear 5.2644e-09 .
                    Dirichet max_error_linear 1.3701e-10, max_error_non_linear 8.628307e-10  max_compare 7.1054e-15.
                    mix_1 , max_error_linear 1.9176e-07, max_error_non_linear 2.3213e-08  max_compare 2.2865e-12.
                    mix_2: max_error_linear 5.3148e-10, max_error_non_linear 7.7868e-10  max_compare 1.5737e-14.
                    %}
                else
                    res_inv = 1;
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
                pos_this = find(x_plot<=e & x_plot>=s);
                x_plot = x_plot(pos_this);
                v_plot = X_inte(2,:);
                u_plot = X_inte(3,:);
                %z_plot = X_inte(4,:);
                v_plot = v_plot(pos_this);
                u_plot = u_plot(pos_this);
                %z_plot = z_plot(pos_this);
                v_true = fun_sol(x_plot);
                %u_true = fun_der_sol(x_plot);
                %z_true = fun_der2_sol(x_plot);
                v_error_inte = v_plot - v_true;
                %u_error_inte = u_plot - u_true;
                %z_error_inte = z_plot - z_true;
                v_error_rk4 = X_rk4(1,:) - v_true;
                %u_error_rk4 = X_rk4(2,:) - u_true;
                %z_error_rk4 = fun([x_plot; X_rk4]) - z_true;
                alpha_inte = A1(1)*v_plot(1)+A1(2)*u_plot(1)+A1(3)*v_plot(end)+A1(4)*u_plot(end);
                beta_inte = A2(2)*v_plot(1)+A2(2)*u_plot(1)+A2(3)*v_plot(end)+A2(4)*u_plot(end);
                alpha_rk4 = A1(1)*X_rk4(1,1)+A1(2)*X_rk4(2,1)+A1(3)*X_rk4(1,end)+A1(4)*X_rk4(2,end);
                beta_rk4 = A2(2)*X_rk4(1,1)+A2(2)*X_rk4(2,1)+A2(3)*X_rk4(1,end)+A2(4)*X_rk4(2,end);
                output_max_error(ii,:) = [ii, max(abs(v_error_inte)), max(abs(v_error_rk4)), max_error_linear, alpha_inte-alpha, alpha_rk4-alpha, beta_inte-beta, beta_rk4-beta,];
                %output_max_error(i,:) = [i, max(abs(v_error_inte)), max(abs(v_error_rk4)), max(abs(u_error_inte)), max(abs(u_error_rk4)), max(abs(z_error_inte)), max(abs(z_error_rk4))];
                %error_change_inte = v_error_inte(2:end)-v_error_inte(1:end-1);
                %error_change_rk4 = v_error_rk4(2:end)-v_error_rk4(1:end-1);
                
                
                
                if abs(myobjopt_inte)< 10^-6
                    %{
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
                    %}
                    output_coef(1:2, ii) = [a_1, a_0];
                    output_coef(3:end, ii) = -bopi^2*J_N(2:M).^2.*b_M(2:M);
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
                z_right_plot = fun(X);
                %z_right_plot = -2*theta*sin(theta*x_R_plot)-theta^2*x_R_plot.*cos(theta*x_R_plot)-(x_R_plot.*cos(theta*x_R_plot)).^2 - 0.1*(cos(theta*x_R_plot)-theta*x_R_plot.*sin(theta*x_R_plot)).^2 + v_M_plot.^2 + 0.1*u_M_plot.^2;
                diff = z_left_plot-z_right_plot;
                pos_this_ = find(x_R_plot<=e & x_R_plot>=s);
                diff = diff(pos_this_);
                if ii == plot_scen 
                    %{
                    filename_fig = strcat('output/fourier_normalizer_app_ode_order2_en_unique_scen_', num2str(i), '_' ,cond_type{1},  '_ps_' , para_scen , '_condi_D_', apply_condition_diri, '_q_', num2str(q), '_deg_', num2str(deg), '_theta_', num2str(round(theta)), '_p_vv_', num2str(p_vv) , '_q_plot_', num2str(q_plot),  figtype_fig);
                    title = strcat('The plot of $(x,y^{(2)}-f(x,y,y^{(1)}))$ over $[s,e]$');
                    legend_y = 'diff';
                    xlabel_this = '$x$';
                    ylabel_this = '';
                    plot_latex(filename_fig, x_R_plot(pos_this),diff, title, legend_y, xlabel_this, ylabel_this, plot_location, 'r')
                   %}
                    v_true_plot = fun_sol(x_R_plot);
                    xlabel_this = '$x$';   
                    filename_fig = strcat('output/fourier_normalizer_app_ode_order2_en_test_yv_scen_', num2str(i), '_' ,cond_type{1}, 'init_neum_', num2str(init_by_neumann) , '_ps_' , para_scen , '_condi_D_', apply_condition_diri, '_q_', num2str(q), '_deg_', num2str(deg), '_theta_', num2str(round(theta)), '_p_vv_', num2str(p_vv), '_q_plot_', num2str(q_plot), figtype_fig);
                    title = strcat('The $v$ and $y$ over $[0,b]$');
                    legend_y = '$v$';
                    legend_z = '$y$';
                    ylabel_this = '';
                    plot_latex_2(filename_fig,x_R_plot, v_M_plot,  v_true_plot, title, legend_y, legend_z,  xlabel_this, ylabel_this, plot_location);
                    filename_true_opt_rk4_fig = strcat('output/fourier_normalizer_app_ode_order2_en_test_true_opt_rk4_scen_', num2str(i), '_' ,cond_type{1}, 'init_neum_', num2str(init_by_neumann) , '_ps_' , para_scen , '_condi_D_', apply_condition_diri, '_q_', num2str(q), '_deg_', num2str(deg), '_theta_', num2str(round(theta)), '_p_vv_', num2str(p_vv), '_q_plot_', num2str(q_plot), figtype_fig);
                    title = strcat('$opt$, $rk4$, $y$ over $[s,e]$');
                    legend_y = '$opt$';
                    legend_w = '$y$';
                    legend_w = '$rk4$';
                    
                    %plot_latex_3(filename_true_opt_rk4_fig, x_plot, v_plot,v_true, X_rk4(1,:), title, legend_y, legend_z, legend_w, xlabel_this, ylabel_this, plot_location )
                    
                    if plot_only
                        return
                    end
                end
            %catch
            %    continue
            %end
            %return
            output(ii,:) = [ii, myobjopt_inte, myobjopt_shooting, max(abs(diff)), v_plot(1)-z_1_true,  u_plot(1)-z_2_true,  v_plot(end)-z_3_true, u_plot(end)-z_4_true,  X_rk4(1,1)-z_1_true, X_rk4(2,1)-z_3_true, X_rk4(1,end)-z_3_true, X_rk4(2,end)-z_4_true];
        end
        vnames = {'scen_id', 'err_opt_inte', 'err_opt_shooting', 'max_diff', 'err_v_s_inte', 'err_u_s_inte', 'err_v_e_inte', 'err_u_e_inte', 'err_v_s_rk4', 'err_u_s_rk4', 'err_v_e_rk4', 'err_u_e_rk4'};
        fourier_normalizer_app_ode_deg_2_uniquenss_file = ['output/fourier_normalizer_app_ode_order2_en_test_cutoff_u_', num2str(cut_off_u), '_cutoff_v_',num2str(cut_off_v), '_', cond_type{1}, 'init_neum_', num2str(init_by_neumann) , '_ps_' , para_scen , '_condi_D_', apply_condition_diri, '_q_', num2str(q), '_deg_', num2str(deg), '_theta_', num2str(round(theta)) , '_p_vv_', num2str(p_vv) , '_q_plot_', num2str(q_plot), '.xlsx'];
        mylib_writearray(vnames, output, fourier_normalizer_app_ode_deg_2_uniquenss_file);
        fourier_normalizer_app_ode_deg_2_max_error_file = ['output/fourier_normalizer_app_ode_order2_en_test_max_error_cutoff_u_', num2str(cut_off_u), '_cutoff_v_',num2str(cut_off_v), '_', cond_type{1}, 'init_neum_', num2str(init_by_neumann) , '_ps_' , para_scen , '_condi_D_', apply_condition_diri, '_q_', num2str(q), '_deg_', num2str(deg), '_theta_', num2str(round(theta)), '_p_vv_', num2str(p_vv), '_q_plot_', num2str(q_plot),'.xlsx'];
        vnames = {'scen_id', 'v_err_inte', 'v_err_rk4', 'linear_error' ,'u_err_inte', 'u_err_rk4', 'z_err_inte', 'z_err_rk4'};
        mylib_writearray(vnames, output_max_error, fourier_normalizer_app_ode_deg_2_max_error_file);
        vnames = {};
        for iii = 1:length(init_guesses_der)
            vnames = [vnames, strcat('scen_',num2str(iii))];
        end
        fourier_normalizer_app_ode_deg_2_coef_file = ['output/fourier_normalizer_app_ode_order2_en_test_coef_cutoff_u_', num2str(cut_off_u), '_cutoff_v_',num2str(cut_off_v), '_', cond_type{1}, 'init_neum_', num2str(init_by_neumann) , '_ps_' , para_scen , '_condi_D_', apply_condition_diri, '_q_', num2str(q), '_deg_', num2str(deg), '_theta_', num2str(round(theta)), '_p_vv_', num2str(p_vv), '_q_plot_', num2str(q_plot),'.xlsx'];
        mylib_writearray(vnames, output_coef, fourier_normalizer_app_ode_deg_2_coef_file);
    end
end

%function V = investigate(X_inte, a_0, a_1, b, G, P, Q, R, alpha, beta, m, n, cond_type,  A1,A2)
function V = investigate(X_inte, a_0, a_1, b, G, Q, R, alpha, beta, m, n, cond_type,  A1,A2)
    %P = [P,0];
    Q = [Q,0];
    R = [R,0];
    z_M = X_inte(4,:);
    M = length(z_M);
    N = 2*M;
    z_N = [0,-fliplr(z_M(2:M)), z_M];  
    J_M = zeros(1,M);
    J_M(2:M) = 1./(1:M-1);
    A = ones(1,N);
    A(2:2:end) = -1;
    b_M = 2*imag(ifft(z_N));
    b_M = b_M.*A;
    b_M = b_M(1:M);
    b_0 = b_M(1);
    b_M = b_M(2:end);
    x = X_inte(1,:);    
    x_R_plusb = zeros(1, M+1);
    x_R_plusb(1:M) = x;
    x_R_plusb(M+1) = b;
    x_0 = x(1);
    x = x(2:end);
    v = X_inte(2,:);
    v_Mplus1 = zeros(1,M+1);
    u_Mplus1 = zeros(1,M+1);
    v_Mplus1(1:M) = v;
    
    v_0 = v(1);
    v_M = a_1 + a_0*b;
    v_Mplus1(M+1) = v_M;
    v = v(2:end);
    u = X_inte(3,:);
    u_Mplus1(1:M) = u;
    u_0 = u(1);
    u = u(2:end);
    z = X_inte(4,:);
    z_0 = z(1);
    z = z(2:end);
    bopi = b/pi;
    A = zeros(M-1,M-1);
    C = zeros(M-1,M-1);
    K = 1:(M-1);
    I = ones(M-1,1);
    alt = I';
    alt(1:2:end)=-1;
    alt = -alt;
    
    %v_M = a_1+b*a_0;
    u_M = a_0 - bopi*(-1)*sum(alt.*b_M.*(1./K));
    u_Mplus1(end) = u_M;
    for i = 1:M-1
        A(i,:) = sin(2*pi*K*i/N);
        C(i,:) = cos(2*pi*K*i/N);
    end
    
    lambda = 2*b/N;
    b_check = 2*a_1/(bopi^2*M)*diag(K.^2)*A*I + 2*lambda*a_0/(M*bopi^2)*diag(K.^2)*A*K' - 2/(M*bopi^2)*diag(K.^2)*A*v';
    b_check_result = b_M'-b_check;
    u_check = v_0*(2/(M^2*bopi)*C*diag(K)*A*K' - I/b-2/(M*bopi)*C*diag(K)*A*I) + v_M*(I/b - 2/(M^2*bopi) *C *diag(K)*A*K') + 2/(bopi*M)*C*diag(K)*A*v';
    u_check_result = u'-u_check;
    

    u_0_check_1 =v_0*(-1/b + 1/(M*bopi)*sum(alt.*K.*cot(pi*K/N)) - 1/bopi*sum(alt.*cot(pi*K/N)));
    u_0_check_2 =  v_M *(1/b - 1/(M*bopi)*sum(alt.*K.*cot(pi*K/N)));
    u_0_check_3 =  1/(bopi)*sum(alt.*v.*cot(K*pi/N));
    

     u_check_0_2_1 = v_0*(2/(M^2*bopi)*I'*diag(K)*A*K' - 1/b-2/(M*bopi)*I'*diag(K)*A*I); 
     u_check_0_2_2 = v_M*(1/b - 2/(M^2*bopi) *I' *diag(K)*A*K');
     u_check_0_2_3 =  2/(bopi*M)*I'*diag(K)*A*v';
    r1 = u_0_check_1 - u_check_0_2_1;
    r2 = u_0_check_2 - u_check_0_2_2;
    r2 = u_0_check_3 - u_check_0_2_3;
    u_0_lim = (v_0-v_M)/(1-2/b);
    
    %v_M = a_1+b*a_0;
    %u_M = a_0 - bopi*(-1)*sum(alt.*b_M.*(1./K));
    
    u_check_M_2_1 = v_0*(2/(M^2*bopi)*(-alt)*diag(K)*A*K' - 1/b-2/(M*bopi)*(-alt)*diag(K)*A*I); 
    u_check_M_2_2 = v_M*(1/b - 2/(M^2*bopi) *(-alt) *diag(K)*A*K');
    u_check_M_2_3 =  2/(bopi*M)*(-alt)*diag(K)*A*v';
    u_M_r = u_check_M_2_1 + u_check_M_2_2 + u_check_M_2_3;
    
    check_u_M_1 = -v_0*(1/b+1/(M*bopi)*sum(K.*cot(K*pi/N))) +v_0*1/(bopi)*sum((alt).*tan(K*pi/N));
    check_u_M_2 = v_M*(1/b+1/(M*bopi)*sum(K.*cot(K*pi/N)));
    check_u_M_3 = -1/(bopi)*sum((alt).*v.*tan(K*pi/N));
    check_u_M_r = check_u_M_1 + check_u_M_2 + check_u_M_3;
    
    u_check_details = zeros(M-1,1);
    for i = 1:M-1
        cot_i = cot((K+i)*pi/N) + cot((K-i)*pi/N);
        cot_i(i) = cot((i+i)*pi/N);
         temp1 = v_0 *(0.5/bopi*sum((-1)^(i+1)*alt.*cot_i) - 0.5/(bopi*M)*sum((-1)^(i+1)*alt.*cot_i.*K) -1/b);
         temp2 = v_M*(0.5/(bopi*M)*sum((-1)^(i+1)*alt.*cot_i.*K)  +1/b );
         temp3 = -0.5/bopi*sum((-1)^(i+1)*alt.*cot_i.*v);
         u_check_details(i,1) = temp1 +temp2 +temp3;
    end
    disp(max(abs(u_check - u_check_details)));
    % U = Alpha*V,  
    Alpha = zeros(1+M, 1+M);
    Alpha(1,1)= 1/(M*bopi)*sum(alt.*K.*cot(K*pi/N))-1/b-1/bopi*sum(alt.*cot(K*pi/N));
    
    Alpha(1,M+1) = -1/(M*bopi)*sum(alt.*K.*cot(K*pi/N))+1/b;
    Alpha(1,2:M) = 1/bopi*alt.*cot(K*pi/N);
    Alpha(M+1,1) = -1/(M*bopi)*sum(K.*cot(K*pi/N))-1/b+1/bopi*sum(alt.*tan(K*pi/N));
    Alpha(M+1,M+1) = 1/b+1/(M*bopi)*sum(K.*cot(K*pi/N));
    Alpha(M+1,2:M) = -1/bopi*alt.*tan(K*pi/N);
    for i = 1:M-1
        cot_i = cot((K+i)*pi/N) + cot((K-i)*pi/N);
        cot_i(i) = cot((i+i)*pi/N);
        Alpha(i+1,1) = 0.5/bopi*sum((-1)^(i+1)*alt.*cot_i) - 0.5/(bopi*M)*sum((-1)^(i+1)*alt.*cot_i.*K) -1/b;
        Alpha(i+1,M+1) = 0.5/(bopi*M)*sum((-1)^(i+1)*alt.*cot_i.*K)  +1/b ; 
        Alpha(i+1,2:M) = -0.5/bopi*(-1)^(i+1)*alt.*cot_i;
    end
    disp(max(abs(u_Mplus1' - Alpha*v_Mplus1')));
    
    O = (2/M)^0.5*A;
    Theta = O*diag(1./(K.^2))*O;

    
    g = zeros(1,M+1);
    g(1) = -alpha;
    g(M+1) = -beta;
    g(2:M) = Theta*G(2:M)';
    g = -g;
    Phi = zeros(M+1,M+1);
    if strcmp(cond_type, 'Dirichlet')
        Phi(1,:) = 0;
        Phi(1,m+1) = 1;
        Phi(1+M,:) = 0;
        Phi(1+M,m+n+1) = 1;
    else
        Phi(1,:) = A1(2)*Alpha(m+1,:) + A1(4)*Alpha(m+n+1,:);
        Phi(1,m+1) = Phi(1,m+1) +  A1(1);
        Phi(1,m+n+1) = Phi(1,m+n+1) +  A1(3);
        Phi(M+1,:) = A2(2)*Alpha(m+1,:) + A2(4)*Alpha(m+n+1,:);
        Phi(M+1,m+1) = Phi(M+1,m+1) +  A2(1);
        Phi(M+1,m+n+1) = Phi(M+1,m+n+1) +  A2(3);
        % first row assolucat to condition 1
    end
    for i = 1:(M-1)
        Phi(i+1,1) = -(M-i)/(M*bopi^2) + (Theta(i,:).*R(2:M))*Alpha(2:M,1);
        Phi(i+1,1+M) = -i/(M*bopi^2) + (Theta(i,:).*R(2:M))*Alpha(2:M,1+M);
        Phi(i+1,2:M) = Theta(i,:).*Q(2:M) +  (Theta(i,:).*R(2:M))*Alpha(2:M,2:M);
        Phi(i+1,i+1) = Phi(i+1,i+1) + 1/bopi^2;
    end
    %V = inv(Phi)*g';
    V = Phi\g';
    %disp(max(abs(V-v_Mplus1')));
    
    %pos = find(x_R_plusb>=1 & x_R_plusb<=3);
    %x_plot_temp = x_R_plusb(pos);
    %v_plot_temp_true = v_R_plusb_true(pos);
    %V_plot_temp = V(pos);
    %v_by_opt = v_Mplus1(pos);
    %max(abs(v_by_opt-v_plot_temp_true));
    %val = max(abs(V-v_Mplus1'));
    %{
    xlabel_this = '$x$';   
    filename_fig = strcat('output/fourier_normalizer_app_ode_order2_en_test_investigate_', figtype_fig);
    title = strcat('The $err_{opt}$ and $err_{linear}$ over $[s,e]$');
    legend_y = '$opt$';
    legend_z = '$linear$';
    ylabel_this = '';
    plot_latex_2(filename_fig,x_plot_temp, v_by_opt-v_plot_temp_true,  V_plot_temp'-v_plot_temp_true, title, legend_y, legend_z,  xlabel_this, ylabel_this, plot_location);
    %}
    
end
    

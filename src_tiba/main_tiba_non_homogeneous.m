
%{
%}

function main_tiba_non_homogeneous()
    global cuf_off_para 
    %cut_off_u = true;
    %cut_off_v = true;
    cut_off_u = false;
    cut_off_v = false;
    flag_display = false;
    %flag_plot = true;
    init_by_neumann = false;
    plot_only = false;
    save_coef = false;
    
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
    q = 7; %default setting in paper
    %q = 8;
    p = q-1; %p should be smaller than q,  
    q_plot = 10; %default used in paper
    
    %deg = 1; %donot change it.  deg=2 might not work
    deg = 2; %donot change it.  deg=2 might not work
    para = 1.5; %or para =1
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
    

    para_scen = 'uv';  %test linear ODE    
    %cond_types = {'Neumann','Dirichlet', 'mix_1', 'mix_2'};
    cond_types = {'Neumann'};
    
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
    

    p_para = {0.0, 0.0, 0, 0.1, 1}; % default  

    [p_uu, p_uv, p_vv, p_u, p_v] = p_para{:};
    function val = fun_sol(X)    %only for x<=0 old:   change to x>=0
        if strcmp(task, 'non_linear')
            x = abs(X);
            if deg==0
                val = cos(theta*x);
            else
                val = x.^deg.*cos(theta*x);
            end
        else
            disp(['not implemented with task ', task])
            return;
        end
    end

    function val = fun_der_sol(X)  
        if strcmp(task, 'non_linear')
            x = abs(X);
            if deg ==0
                val = -theta*sin(theta*x);
            elseif deg==1
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
            if deg == 0
                val = -theta^2*cos(theta*x);
            elseif deg==1
                val = -2*theta*sin(theta*x)-theta^2*x.*cos(theta*x);
            elseif deg==2
                val = deg*(deg-1)*cos(theta*x)-2*deg*theta*x.*sin(theta*x)-theta^2*x.^2.*cos(theta*x);
            else
                val = deg*(deg-1)*x.^(deg-2).*cos(theta*x)-2*deg*theta*x.^(deg-1).*sin(theta*x)-theta^2*x.^deg.*cos(theta*x);
            end
        else
            disp(['not implemented with task ', task])
            return;
        end
    end

    function val = fun(X)
        if strcmp(task, 'non_linear')
            x = X(1,:);
            v = X(2,:);
            u = X(3,:);
            if deg == 0
                der_2 = -theta^2*cos(theta*x);
            elseif deg==1
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
        else
            disp(['not implemented with task ', task])
            return;
        end
    end
    
    function val = generate_r(x_R)
        z =  fun_der2_sol(x_R);
        v = fun_sol(x_R);
        u = fun_der_sol(x_R);
        %val = z-p_vv*v.^2-p_v*v-p_u*u;
        val = z-p_v*v-p_u*u; %only linear term!
    end

    function val = fun_q(X)
        %val = p_v*ones(1,length(X)).*co_R;
        val = p_v*ones(1,length(X));
    end

    function val = fun_p(X)
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
        %randoms = rand(2,5);
        randoms = rand(2,1);
        init_guesses_der = randoms(1,:)-0.5;
        init_guesses_fun = randoms(2,:)-0.5;
        init_guesses_der = [cheat_derivative + init_guesses_der, 2*cheat_derivative + init_guesses_der, -2*cheat_derivative + init_guesses_der, 3*cheat_derivative + init_guesses_der, -3*cheat_derivative + init_guesses_der];
        init_guesses_fun = [cheat_fun + init_guesses_fun, 2*cheat_fun + init_guesses_fun, -2*cheat_fun + init_guesses_fun, 5*cheat_fun + init_guesses_fun, -5*cheat_fun + init_guesses_fun];

        output = zeros(length(init_guesses_der), 9);
        output_max_error = zeros(length(init_guesses_der),11);
        output_coef = zeros(M+1, length(init_guesses_der));
        if strcmp(cond_type, 'Dirichlet')
            A1 = [1,0,0,0];
            A2 = [0,0,1,0]; 
        elseif strcmp(cond_type, 'Neumann')
            A1 = [1,0,0,0];
            A2 = [0,1,0,0]; 
        elseif strcmp(cond_type, 'mix_1')
            A1 = [1,0,0,0];
            A2 = [0,0,0,1];
        elseif strcmp(cond_type, 'mix_2')
            A1 = [1,1,0,0];
            A2 = [0,0,1,1];
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

            %x_M = X_inte(1,:);
            Q = fun_q(x_R);
            Q = Q.*co_R;
            P = fun_p(x_R);
            P = P.*co_R;
            R = generate_r(x_R);
            R = R.*co_R;
            max_error_linear = 0;
            if  strcmp(para_scen, 'uv') && abs(myobjopt_inte)<10^(-6)
                V_linear =  engine_tiba(b, R, Q, P, alpha, beta, m, n, cond_type,  A1,A2);   
                x_R_plusb = zeros(1, M+1);
                x_R_plusb(1:M) = x_R;
                x_R_plusb(M+1) = b;
                v_true_this = fun_sol(x_R_plusb);
                pos_this_ = find(x_R_plusb>=s & x_R_plusb<=e);
                max_error_linear = max(abs(v_true_this(pos_this_)'-V_linear(pos_this_)));                   
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
            v_plot = v_plot(pos_this);
            u_plot = u_plot(pos_this);
            v_true = fun_sol(x_plot);
            v_error_inte = v_plot - v_true;
            v_error_rk4 = X_rk4(1,:) - v_true;
            alpha_inte = A1(1)*v_plot(1)+A1(2)*u_plot(1)+A1(3)*v_plot(end)+A1(4)*u_plot(end);
            beta_inte = A2(1)*v_plot(1)+A2(2)*u_plot(1)+A2(3)*v_plot(end)+A2(4)*u_plot(end);
            alpha_rk4 = A1(1)*X_rk4(1,1)+A1(2)*X_rk4(2,1)+A1(3)*X_rk4(1,end)+A1(4)*X_rk4(2,end);
            beta_rk4 = A2(1)*X_rk4(1,1)+A2(2)*X_rk4(2,1)+A2(3)*X_rk4(1,end)+A2(4)*X_rk4(2,end);
            if abs(myobjopt_inte)< 10^-6
                output_coef(1:2, ii) = [a_1, a_0];
                output_coef(3:end, ii) = -bopi^2*J_N(2:M).^2.*b_M(2:M);
            end
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
            diff = z_left_plot-z_right_plot;
            pos_this_ = find(x_R_plot<=e & x_R_plot>=s);
            diff = diff(pos_this_);
            output_max_error(ii,:) = [ii, max(abs(v_error_inte)), max(abs(v_error_rk4)), max_error_linear, max(abs(diff)), myobjopt_inte, myobjopt_shooting, alpha_inte-alpha, alpha_rk4-alpha, beta_inte-beta, beta_rk4-beta];
            if ii == plot_scen 
                v_true_plot = fun_sol(x_R_plot);
                xlabel_this = '$x$';   
                filename_fig = strcat('output/tiba_non_homo_yv_scen_', num2str(i), '_' ,cond_type{1}, 'init_neum_', num2str(init_by_neumann) , '_ps_' , para_scen , '_condi_D_', apply_condition_diri, '_q_', num2str(q), '_deg_', num2str(deg), '_theta_', num2str(round(theta)), '_p_vv_', num2str(p_vv), '_q_plot_', num2str(q_plot), figtype_fig);
                title = strcat('The $v$ and $y$ over $[0,b]$');
                legend_y = '$v$';
                legend_z = '$y$';
                ylabel_this = '';
                plot_latex_2(filename_fig,x_R_plot, v_M_plot,  v_true_plot, title, legend_y, legend_z,  xlabel_this, ylabel_this, plot_location);                    
                if plot_only
                    return
                end
            end
            output(ii,:) = [ii, v_plot(1)-z_1_true,  u_plot(1)-z_2_true,  v_plot(end)-z_3_true, u_plot(end)-z_4_true,  X_rk4(1,1)-z_1_true, X_rk4(2,1)-z_3_true, X_rk4(1,end)-z_3_true, X_rk4(2,end)-z_4_true];
        end
        vnames = {'scen_id', 'tibo_err_v_s', 'tibo_err_u_s', 'tibo_err_v_e', 'tibo_err_u_e', 'rk4_err_v_s', 'rk4_err_u_s', 'rk_err_v_e', 'rk4_err_u_e'};
        fourier_normalizer_app_ode_deg_2_uniquenss_file = ['output/tiba_non_homo_cutoff_u_', num2str(cut_off_u), '_cutoff_v_',num2str(cut_off_v), '_', cond_type{1}, 'init_neum_', num2str(init_by_neumann) , '_ps_' , para_scen , '_condi_D_', apply_condition_diri, '_q_', num2str(q), '_deg_', num2str(deg), '_theta_', num2str(round(theta)) , '_p_vv_', num2str(p_vv) , '_q_plot_', num2str(q_plot), '.xlsx'];
        mylib_writearray(vnames, output, fourier_normalizer_app_ode_deg_2_uniquenss_file);
        fourier_normalizer_app_ode_deg_2_max_error_file = ['output/tiba_non_homo_max_error_cutoff_u_', num2str(cut_off_u), '_cutoff_v_',num2str(cut_off_v), '_', cond_type{1}, 'init_neum_', num2str(init_by_neumann) , '_ps_' , para_scen , '_condi_D_', apply_condition_diri, '_q_', num2str(q), '_deg_', num2str(deg), '_theta_', num2str(round(theta)), '_p_vv_', num2str(p_vv), '_q_plot_', num2str(q_plot),'.xlsx'];
        vnames = {'scen_id', 'tibo_error', 'rk4_err', 'linear_err' , 'max_diff', 'fval', 'fval_shooting', 'alpha_err_tibo', 'alpha_err_rk4', 'beta_err_tibo', 'beta_err_rk4'};
        mylib_writearray(vnames, output_max_error, fourier_normalizer_app_ode_deg_2_max_error_file);

        if save_coef
            vnames = {};
            for iii = 1:length(init_guesses_der)
                vnames = [vnames, strcat('scen_',num2str(iii))];
            end
            fourier_normalizer_app_ode_deg_2_coef_file = ['output/tiba_non_homo_coef_cutoff_u_', num2str(cut_off_u), '_cutoff_v_',num2str(cut_off_v), '_', cond_type{1}, 'init_neum_', num2str(init_by_neumann) , '_ps_' , para_scen , '_condi_D_', apply_condition_diri, '_q_', num2str(q), '_deg_', num2str(deg), '_theta_', num2str(round(theta)), '_p_vv_', num2str(p_vv), '_q_plot_', num2str(q_plot),'.xlsx'];
            mylib_writearray(vnames, output_coef, fourier_normalizer_app_ode_deg_2_coef_file);
        end
    end
end
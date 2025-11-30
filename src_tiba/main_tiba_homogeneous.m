
%{
test driver for the following paper

"Trigonometric Interpolation Based Based Approach for Second Order ODE with
Mixed Boundary Conditions"

how to run
keep the default setting except

test_no_sol = false;  %run the default result with solutions,  either
unqiue solution or infinite many solution

test_no_sol = true;  %run the results without solution

%}

function main_tiba_homogeneous()
    global cuf_off_para 
    cut_off_u = false;
    cut_off_v = false;
    flag_display = false;
    test_no_sol = false;
    init_by_neumann = false;
    %flag_plot = true;
    save_coef = false;
    save_bdy = false;
    
    %apply_condition_diri = 'upbound';
    apply_condition_diri = 'no';
    %apply_condition_diri = 'lowbound'; %'no', 'upbound'
    apply_condition_on_lower_boundary_flag = false;
    apply_condition_on_lower_boundary_value = -0.1;
    apply_condition_on_upper_boundary_flag = false;
    apply_condition_on_upper_boundary_value = 10^6;
    
    
    figtype_fig = '.fig';
    %task = 'linear';
    %task = 'non_linear';
    %init_method = 'runge_kutta';
    plot_location = 'northwest';
    constrainttolerance = 10^-15;
    maxiter = 3000;

    fourier_normalizer_const()
    s = 1;
    e = 3*s;
    q = 7; %default setting in paper
    p = q-1; %p should be smaller than q,  
    
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
    
    %cond_types = {'Neumann'};
    %cond_types = {'Dirichlet'};
    %cond_types = {'mix_1'};
    %cond_types = {'mix_2'};

    cond_types = {'Neumann', 'Dirichlet', 'mix_1', 'mix_2'};
    
    
    
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
    
    c_1 = 1;
    c_2 = 1;
    p_v = -1.25*pi^2;
    p_u = -2*pi;
    
    
    function val = fun_sol(X)    %only for x<=0 old:   change to x>=0
        x = abs(X);
        temp1 = exp(-(x-1)*pi).*(cos(pi*(x-1)/2) + 2*sin(pi*(x-1)/2));
        temp2 = exp(-(x-1)*pi).*sin(pi*(x-1)/2);
        val = c_1*temp2 + c_2*temp1;
    end

    function val = fun_der_sol(X)  
        x = abs(X);
        temp1_der = -2.5*pi*exp(-(x-1)*pi).*sin(pi*(x-1)/2);
        temp2_der = exp(-(x-1)*pi).*(0.5*pi*cos(pi*(x-1)/2) -pi*sin(pi*(x-1)/2));
        val = c_1*temp2_der + c_2*temp1_der;
    end

    %F(x,y',y):F(x,u,v)
    %F(x,u,v) = -2*pi*u - 1.25*pi^2*v
    function val = fun(X)
        %x = X(1,:);
        v = X(2,:);
        u = X(3,:);
        val = p_u*u + p_v*v;
    end

    function val = fun_grid(X)   
        val = fun(X);
        if length(val)==length(co_R)
            val = val.*co_R;
        else
            x1 = X(1,:);
            t=length(val);
            for i=1:t
                pos_this = find(x_R>=x1(i));
                val(i) =val(i)*co_R(pos_this(1));
            end
        end
    end

    function val = partial_u(X) %X=(x,v,u)
        u = X(3,:);
        if cut_off_u
            c_u = fourier_normalizer_cut_off(u, l_u,s_u,e_u,r_u,cuf_off_para);
            c_d_u = fourier_normalizer_cut_off_der1(u, l_u,s_u,e_u,r_u,cuf_off_para);
        else
            c_u = ones(1, length(u));
            c_d_u = zeros(1, length(u));
        end
        %val = 0.2*u;
        val =  p_u*(c_u + u.*c_d_u);
    end

    function val = partial_v(X) %X=(x,v,u)
        v = X(2,:);
        if cut_off_v
            c_v = fourier_normalizer_cut_off(v, l_v,s_v,e_v,r_v,cuf_off_para);
            c_d_v = fourier_normalizer_cut_off_der1(v, l_v,s_v,e_v,r_v,cuf_off_para);
        else
            c_v = ones(1, length(v));
            c_d_v = zeros(1, length(v));
        end
        val =  p_v*(c_v + v.*c_d_v);
    end
    % z = g(x) + q*v*h(x) + r*u*h(x)

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
    cond_b = 0;

    for jj = 1:length(cond_types)
        cond_type = cond_types(jj);
        rng('default')
        randoms = rand(2,1);
        cheat_derivative = z_2_true;
        cheat_fun = z_1_true;
        init_guesses_der = randoms(1,:)-0.5;
        init_guesses_fun = randoms(2,:)-0.5;
        
        init_guesses_der = [cheat_derivative + init_guesses_der, 2*cheat_derivative + init_guesses_der, -2*cheat_derivative + init_guesses_der, 3*cheat_derivative + init_guesses_der, -3*cheat_derivative + init_guesses_der];
        init_guesses_fun = [cheat_fun + init_guesses_fun, 2*cheat_fun + init_guesses_fun, -2*cheat_fun + init_guesses_fun, 5*cheat_fun + init_guesses_fun, -5*cheat_fun + init_guesses_fun];
        
        if strcmp(cond_type, 'Neumann')
            init_guesses_der = [cheat_derivative];
            init_guesses_fun = [init_guesses_fun];
        end
        
        %init_guesses_der = [cheat_derivative + init_guesses_der];
        %init_guesses_fun = [cheat_fun + init_guesses_fun];
        output = zeros(length(init_guesses_der), 9);        
        output_max_error = zeros(length(init_guesses_der), 11);
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
            A1 = [1,0,0,0];
            A2 = [0,0,0,1];
            %alpha = fun_sol(s);
            %beta = fun_der_sol(e);
        elseif strcmp(cond_type, 'mix_2')
            A1 = [1,1,0,0];
            A2 = [0,0,1,1];
        end

        alpha = A1(1)*z_1_true+A1(2)*z_2_true+A1(3)*z_3_true+A1(4)*z_4_true;
        beta = A2(1)*z_1_true+A2(2)*z_2_true+A2(3)*z_3_true+A2(4)*z_4_true;
        if test_no_sol
            beta = beta*1.1; %force y_b is not a solution,  so the strutrue between alpha and beta by y_b failed. 
        end
        AA = [A1;A2];
        greeks = [alpha, beta];
        for ii = 1:length(init_guesses_der)
            if ii~=11
                %continue
            end
            init_guess = [init_guesses_fun(ii),init_guesses_der(ii)];             
            [init_opt, myobjopt_shooting, X_rk4] = fourier_normalizer_app_ode_get_init_value_by_rk4_order2(@fun, s, e,  greeks, AA, n, init_guess, maxiter, constrainttolerance); 
            [~, init_z] = runge_kutta_ode_order_k(@fun_grid, x_R,init_opt, m+1, []); 
            if init_by_neumann
                AA_n  = [1,0,0,0;0,1,0,0];                 
                [X_inte_n,~,~,~,~,~] = fourier_normalizer_app_ode_order2_en(@fun, @partial_u, @partial_v, s, e, p,q, init_opt, AA_n,...
                    flag_display,  constrainttolerance,maxiter, apply_condition_diri,cond_b,...
                    apply_condition_on_lower_boundary_flag, apply_condition_on_lower_boundary_value, apply_condition_on_upper_boundary_flag, apply_condition_on_upper_boundary_value, init_z);
                init_z = X_inte_n(4,:);
            end
            [X_inte, a_0, a_1, b,o,myobjopt_inte] = fourier_normalizer_app_ode_order2_en(@fun, @partial_u, @partial_v, s, e, p,q, greeks, AA,...
                flag_display,  constrainttolerance,maxiter, apply_condition_diri,cond_b,...
                apply_condition_on_lower_boundary_flag, apply_condition_on_lower_boundary_value, apply_condition_on_upper_boundary_flag, apply_condition_on_upper_boundary_value, init_z);

            x_M = X_inte(1,:);
            Q = fun_q(x_M);
            Q = Q.*co_R;  % we should have this F =f*h = h(qy +ry')
            P = fun_p(x_M);
            P = P.*co_R;
            R = zeros(1,length(X_inte(4,:)));
            %V_linear = investigate(X_inte, a_0, a_1, b, G, P, Q, R,  alpha, beta, m, n, cond_type,  A1,A2);
            V_linear = engine_tiba(b, R, Q, P,  alpha, beta, m, n, cond_type,  A1,A2);
            x_this = X_inte(1,:);    
            x_R_plusb = zeros(1, M+1);
            x_R_plusb(1:M) = x_this;
            x_R_plusb(M+1) = b;
            v_true_this = fun_sol(x_R_plusb);
            pos_this_ = find(x_R_plusb>=s & x_R_plusb<=e);
            max_error_linear = max(abs(v_true_this(pos_this_)'-V_linear(pos_this_)));

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
            pos_this_ = find(x_plot<=e & x_plot>=s);
            x_plot = x_plot(pos_this_);
            v_plot = X_inte(2,:);
            u_plot = X_inte(3,:);
            %z_plot = X_inte(4,:);
            v_plot = v_plot(pos_this_);
            u_plot = u_plot(pos_this_);
            %z_plot = z_plot(pos_this);
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


            q_plot = 10;
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

            if ii == 11
                filename_fig = strcat('output/tiba_homo_unique_scen_', num2str(ii), '_' ,cond_type{1}, '_test_no_sol_', num2str(test_no_sol), '_init_neum_', num2str(init_by_neumann) , '_condi_D_', apply_condition_diri, '_q_', num2str(q), '_q_plot_', num2str(q_plot),  figtype_fig);
                title = strcat('The plot of $(x,y^{(2)}-f(x,y,y^{(1)}))$ over $[s,e]$');
                legend_y = 'diff';
                xlabel_this = '$x$';
                ylabel_this = '';
                plot_latex(filename_fig, x_R_plot(pos_this_),diff, title, legend_y, xlabel_this, ylabel_this, plot_location, 'r')
                v_true_plot = fun_sol(x_R_plot);
                xlabel_this = '$x$';   
                filename_fig = strcat('output/tiba_homo_yv_scen_', num2str(ii), '_' ,cond_type{1}, '_tnl_', num2str(test_no_sol), '_q_', num2str(q), figtype_fig);
                title = strcat('The $v$ and $y$ over $[0,b]$');
                legend_y = '$v$';
                legend_z = '$y$';
                ylabel_this = '';
                plot_latex_2(filename_fig,x_R_plot, v_M_plot,  v_true_plot, title, legend_y, legend_z,  xlabel_this, ylabel_this, plot_location);
            end
            output_max_error(ii,:) = [ii, max(abs(v_error_inte)), max(abs(v_error_rk4)), max_error_linear, max(abs(diff)), myobjopt_inte, myobjopt_shooting, alpha_inte-alpha, alpha_rk4-alpha, beta_inte-beta, beta_rk4-beta];
            output(ii,:) = [ii, v_plot(1)-z_1_true,  u_plot(1)-z_2_true,  v_plot(end)-z_3_true, u_plot(end)-z_4_true,  X_rk4(1,1)-z_1_true, X_rk4(2,1)-z_3_true, X_rk4(1,end)-z_3_true, X_rk4(2,end)-z_4_true];
        end
        
        if save_bdy
            vnames = {'scen_id', 'tibo_err_v_s', 'tibo_err_u_s', 'tibo_err_v_e', 'tibo_err_u_e', 'rk4_err_v_s', 'rk4_err_u_s', 'rk_err_v_e', 'rk4_err_u_e'};
            fourier_normalizer_app_ode_deg_2_uniquenss_file = ['output/tiba_homo_', cond_type{1}, '_tns_', num2str(test_no_sol), '_q_', num2str(q),  '_c2_', num2str(c_2), '.xlsx'];
            mylib_writearray(vnames, output, fourier_normalizer_app_ode_deg_2_uniquenss_file);
        end
        fourier_normalizer_app_ode_deg_2_max_error_file = ['output/tiba_homo_max_error_', cond_type{1}, '_test_no_sol_', num2str(test_no_sol), '_q_', num2str(q), '_c2_', num2str(c_2),'.xlsx'];        
        vnames = {'scen_id', 'tibo_error', 'rk4_err', 'linear_err' , 'max_diff', 'fval', 'fval_shooting', 'alpha_err_tibo', 'alpha_err_rk4', 'beta_err_tibo', 'beta_err_rk4'};
        mylib_writearray(vnames, output_max_error, fourier_normalizer_app_ode_deg_2_max_error_file);
        if save_coef
            vnames = {};
            for ii = 1:length(init_guesses_der)
                vnames = [vnames, strcat('scen_',num2str(ii))];
            end
            fourier_normalizer_app_ode_deg_2_coef_file = ['output/tiba_homo_coef_', cond_type{1}, '_tns_', num2str(test_no_sol), '_ps_' , para_scen , '_q_', num2str(q),  '.xlsx'];
            mylib_writearray(vnames, output_coef, fourier_normalizer_app_ode_deg_2_coef_file);
        end
       
    end
end


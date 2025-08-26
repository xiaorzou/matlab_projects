%{
driver for ode_order_1,  optimization with z=u' as variables to be
optimized.
Aug, 3, 2025
%}

function ode_order_1_by_z_main()
    fourier_normalizer_const()
    global cuf_off_para 
    global test_type weight_constrain s  e dim ps qs c3 degs init_method qqs paras ws pps
    global maxfunevals flag_applying_constrains flag_plot figtype_fig plot_location
    
    output = zeros(length(flag_applying_constrains)*length(qqs)*length(paras)*length(degs), 9);
    counter = 1;
    ode_order_1_by_u_enhance_file = strcat('output/ode_order_1_by_z_main_init_', init_method, '_max_fun_', num2str(maxfunevals),'.xlsx');
    vnames = {'constrain', 'q', 'deg', 'theta', 'fval', 'err_tibo', 'err_rk', 'H_tibo', 'H_rk'};
        
        
    for flag_applying_constrain = flag_applying_constrains
        for q = qqs
            p = q-1; %p should be smaller than q,  
            for para = paras
                for deg_x = degs
                    theta = (0.5+para)*pi;
                    if flag_applying_constrain
                        pps(3)=1;  %not work if we set pps(3)=0 due to not sufficient variables.  notice that we use z as variable, not y.  maybe we can consider use y as variable, rather than z.  
                        ws(3) = 1; %use different weights might deterioate model perforamnce.  is it related to some symmetric property 
                    end                    
                    M = 2^q;
                    N = 2*M;
                    n = 2^p; % we should make n as large as possible,  as such, p=q-1 should be always true
                    lambda = (e-s)/n;
                    m = (M-n)/2;
                    delta = lambda*m;
                    o = s - delta;
                    b = e-s + 2*delta;
                    init_pos = m+1;
                    pos_grid = [m+1:m+n+1]; 

                    x_N = -b + o + lambda*(0:N-1);
                    x_M = x_N(1:M);
                    x_R = x_N(M+1:end);
                    co_R = fourier_normalizer_cut_off(x_R, s-delta,s,e,e+delta,cuf_off_para);
                    co_L = [0,fliplr(co_R(2:M))];
                    co_R = co_R';
                    co_der_R = fourier_normalizer_cut_off_der1(x_R, s-delta,s,e,e+delta,cuf_off_para);
                    co_der_R = co_der_R';

                    q_2N = q+1;
                    M_2N = 2^q_2N;
                    N_2N = 2*M_2N;
                    p_2N = q_2N-1;
                    n_2N = 2^p_2N; % we should make n as large as possible,  as such, p=q-1 should be always true
                    lambda_2N = (e-s)/n_2N;
                    x_2N = -b + o + lambda_2N*(0:N_2N-1);
                    x_2N_R = x_2N(M_2N+1:end);
                    co_2N_R = fourier_normalizer_cut_off(x_2N_R, s-delta,s,e,e+delta,cuf_off_para);
                    co_2N_R = co_2N_R';
                    co_2N_der_R = fourier_normalizer_cut_off_der1(x_2N_R, s-delta,s,e,e+delta,cuf_off_para);
                    co_2N_der_R = co_2N_der_R';
                    y_s_vect = fun_sol(s, x_R, test_type, M, co_R, theta, deg_x, c3);                  
                    u_rk = runge_kutta_ode_deg_1_variable_k(@fun_rk, x_R',y_s_vect, init_pos, x_2N_R, co_2N_R, co_2N_der_R, test_type, theta,deg_x, ps, qs, c3);
                    z_rk = fun([x_R', u_rk],x_R, test_type, M, co_R, co_der_R,theta, deg_x, ps, qs, c3);
                    u_rk = u_rk';
                    z_rk = z_rk';
                    u_rk_left = zeros(M, dim);
                    z_rk_left = zeros(M, dim);
                    for i = 1:dim
                        u_rk_left(:,i) =  [0,fliplr(u_rk(i,2:end))]'; 
                        z_rk_left(:,i) =  [0,-fliplr(z_rk(i,2:end))]';
                    end
                    u_rk = u_rk';
                    z_rk = z_rk';
                    H_rk = get_H(u_rk(pos_grid,:),test_type,c3);
                    H_rk_max = max(abs(H_rk));

                    u_real = fun_sol(x_R',x_R, test_type, M, co_R, theta, deg_x, c3);
                    %H_real = get_H(u_real(pos_grid,:),test_type,c3);
                    
                    
                    
                    
   
                    if strcmp(test_type,'base')
                        true_fun_value = fun_sol(x_R', x_R, test_type, M, co_R, theta, deg_x, c3);

                        true_der_value = fun_der_sol(x_R',  x_R,M, co_R, co_der_R, theta, deg_x, c3);
                        true_fun_value = true_fun_value';
                        true_der_value = true_der_value';
                        true_fun_value_left = zeros(M, dim);
                        true_der_value_left = zeros(M, dim);
                        for i = 1:dim
                            true_fun_value_left(:,i) =  [0,fliplr(true_fun_value(i,2:end))]'; 
                            true_der_value_left(:,i) =  [0,-fliplr(true_der_value(i,2:end))]';
                        end
                    end
                    if strcmp(init_method, 'cheat')
                        init_guess = true_der_value_left;

                    elseif strcmp(init_method, 'simple')
                        %u_init = u_true;
                        z_s_vect = fun([s, y_s_vect],x_R, test_type, M, co_R, co_der_R,theta, deg_x, ps, qs, c3); %
                        init_guess = zeros(M,3);
                        u_init = zeros(M,3);
                        index_0 = m+n+1;
                        init_guess(index_0,:) = -z_s_vect;   %z is odd!!!,  m+n+1 relfects to -s
                        u_init(index_0,:) = y_s_vect; 
                        for i = (index_0+1):M
                            for ii =1:3
                                u_init(i, ii) = (init_guess(i-1, ii)*lambda + u_init(i-1, ii))*co_L(i); %u is even!
                            end
                            init_guess(i, :) = -fun([-x_M(i), u_init(i,:)],x_R, test_type, M, co_R, co_der_R,theta, deg_x, ps, qs, c3); % z is odd,  fun is only defined over positive range!!!
                            for ii = 1:3
                                init_guess(i, ii) = init_guess(i, ii)*co_L(i);
                            end
                        end
                        for i = (index_0-1):-1:2
                            for ii = 1:3
                                u_init(i,ii) = (init_guess(i+1, ii)*(-lambda) + u_init(i+1, ii))*co_L(i); %u is even!
                            end
                            init_guess(i,:) = -fun([-x_M(i), u_init(i,:)],x_R, test_type, M, co_R, co_der_R,theta, deg_x, ps, qs, c3)*co_L(i); %z is odd,  fun is only defined over positive range!!!
                            for ii = 1:3
                                init_guess(i,ii) = init_guess(i,ii)*co_L(i); %z is odd!
                            end
                        end

                    elseif strcmp(init_method, 'rk')       
                        init_guess = z_rk_left;
                    else
                        z_s_vect = fun([s, y_s_vect],x_R, test_type, M, co_R, co_der_R,theta, deg_x, ps, qs, c3); %
                        init_guess = zeros(M, 3);
                        for ii = 1:3
                            init_guess(:, ii) = -z_s_vect(ii)*co_L;
                        end
                    end

                    r_val = get_r_value(x_R',x_R, test_type, M, co_R, co_der_R,theta, deg_x, ps, qs, c3);
                    [u_N_out, fval] = ode_order_1_by_z_engine(@fun_w_r, @partial_u,@get_H, @get_H_der1, ...
                        y_s_vect, s, e, p, init_guess, ws, pps, flag_applying_constrain, weight_constrain, r_val, test_type, ps, qs, c3, M, co_R);
                    u_grid = u_N_out(M+1:N,:);
                    u_grid = u_grid(pos_grid,:);
                    
                    u_real_grid = u_real(pos_grid,:);     
                    u_rk_grid = u_rk(pos_grid,:);
                    err_rk_grid = u_rk_grid-u_real_grid;
                    err_tibo_grid = u_grid-u_real_grid;                                                        
                    max_u = max(abs(u_real_grid));
                    for a = 1:dim
                        err_tibo_grid(:,a) = err_tibo_grid(:,a)/max_u(a);
                        err_rk_grid(:,a) = err_rk_grid(:,a)/max_u(a);
                    end
                    err_tibo = max(max(abs(err_tibo_grid)));
                    err_rk_grid = max(max(abs(err_rk_grid)));
                    H_tri = get_H(u_grid, test_type, c3);
                    H_tibo_grid_max = max(abs(H_tri));
                                                    
                    output(counter, :)=[flag_applying_constrain, q, deg_x, theta, fval, err_tibo, err_rk_grid, H_tibo_grid_max, H_rk_max];   
                    counter = counter +1;
                    if flag_plot
                        coef = zeros(M,dim);
                        for a = 1:dim    
                            coef_this = fourier_normalizer_value2coef(u_N_out(:,a)', 'cos');
                            coef_this = coef_this(1:M);
                            coef(:, a)= coef_this;
                        end
                        
                        q_plot = 10;
                        p_plot = q_plot-1; %p should be smaller than q,  
                        M_plot = 2^q_plot;
                        N_plot = 2*M_plot;
                        n_plot = 2^p_plot; 
                        lambda_plot = (e-s)/n_plot;
                        x_N_plot = -b + o + lambda_plot*(0:N_plot-1);
                        x_R_plot = x_N_plot(M_plot+1:end);
                        u_true_plot =  fun_sol_y(x_R_plot, test_type,theta, deg_x, c3);
                        u_tibo_plot = zeros(M_plot, 3);
                        for alpha = 1:3
                            u_R = cos_approx_engine_coef2value(coef(:,alpha), x_R_plot, b);
                            u_tibo_plot(:,alpha) = u_R;
                            u_true_this = u_true_plot(:,alpha);
                            filename_fig = ['output/ode_order_1_by_z_main_alpha_',  num2str(alpha), '_q_', num2str(q), '_t_', num2str(round(theta)), '_deg_', num2str(round(deg_x)), init_method ,  figtype_fig];
                            title = sprintf(char('$\\hat{y}_{%i}$ vs $u_{%i}$  over [0,b]'), alpha, alpha);
                            legend_y = sprintf(char('$\\hat{y}_{%i}$'), alpha);
                            legend_z = sprintf(char('$u_{%i}$'), alpha);
                            xlabel_this = '$x$';
                            ylabel_this = '$y$';
                            plot_latex_2(filename_fig, x_R_plot, u_true_this', u_R, title, legend_y, legend_z, xlabel_this, ylabel_this, plot_location )    
                        end
                    end                    
                end
            end
        end
    end
    if ~flag_plot
        mylib_writearray(vnames, output, ode_order_1_by_u_enhance_file);
    end
end

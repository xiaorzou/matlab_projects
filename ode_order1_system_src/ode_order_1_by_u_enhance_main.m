%{
driver for ode_order_1,  optimization with z=u' as variables to be
optimized.
Aug, 3, 2025
%}

function ode_order_1_by_u_enhance_main()
    fourier_normalizer_const()
    global cuf_off_para 
    global figtype_fig flag_applying_constrains test_type weight_constrain plot_location s  e dim ps qs c3 degs init_method qqs paras flag_plot  ws pps
    output = zeros(length(flag_applying_constrains)*length(qqs)*length(paras)*length(degs), 9);
    counter = 1;
    ode_order_1_by_u_enhance_file = strcat('output/ode_order_1_by_u_enhance_init_', init_method, '.xlsx');
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
                    H_rk = get_H(u_rk(pos_grid,:),test_type,c3);
                    H_rk_max = max(abs(H_rk));

                    u_real = fun_sol(x_R',x_R, test_type, M, co_R, theta, deg_x, c3);
                    H_real = get_H(u_real(pos_grid,:),test_type,c3);

                    u_rk_grid = max(max(abs((u_real(pos_grid)-u_rk(pos_grid)))));

                    true_fun_value = fun_sol(x_R', x_R, test_type, M, co_R, theta, deg_x, c3);

                    true_der_value = fun_der_sol(x_R', M, co_R, co_der_R, theta, deg_x, c3);
                    true_fun_value = true_fun_value';
                    true_der_value = true_der_value';
                    true_fun_value_left = zeros(M, dim);
                    true_der_value_left = zeros(M, dim);
                    for i = 1:dim
                        true_fun_value_left(:,i) =  [0,fliplr(true_fun_value(i,2:end))]'; 
                        true_der_value_left(:,i) =  [0,-fliplr(true_der_value(i,2:end))]';
                    end


                    if strcmp(init_method, 'cheat')
                        init_guess = [true_fun_value_left; zeros(1,dim)];
                        init_guess = [init_guess(1:m+n,:); init_guess(m+n+2:end,:)];

                    elseif strcmp(init_method, 'rk')       
                        init_guess = [u_rk_left; zeros(1,dim)];
                        init_guess = [init_guess(1:m+n,:); init_guess(m+n+2:end,:)];
                    else
                        init_guess = zeros(M, dim);
                        for ii = 1:dim
                            init_guess(:, ii) = y_s_vect(ii)*co_L;
                        end
                    end

                    r_val = get_r_value(x_R',x_R, test_type, M, co_R, co_der_R,theta, deg_x, ps, qs, c3);
                    [u_tibo, fval] = ode_order_1_by_u_enhance_engine(@fun_w_r, @partial_u,@get_H, @get_H_der1, ...
                        y_s_vect, s, e, p, init_guess, ws, pps, flag_applying_constrain, weight_constrain, r_val, test_type, ps, qs, c3, M, co_R);

                    coef_u = zeros(M,dim);
                    alt = ones(1,M);

                    alt(2:2:end)=-1;
                    alt_N = ones(1,N);
                    alt_N(2:2:end)=-1;
                    adjustment = zeros(1,dim);
                    for a = 1:dim
                        u = u_tibo(:, a)';
                        u = [u, fliplr(u(2:M))];     
                        coef_u_a =  2*real(ifft(u));
                        coef_u_a = coef_u_a(1:M).*alt;
                        coef_u_a(1) = sum(u(1:2:end))/M;
                        coef_u(:,a) = coef_u_a;
                        adjustment(a) = sum(u.*alt_N)/M;
                    end

                    u_tibo = fliplr(u_tibo')';
                    u_tibo = u_tibo(1:M,:);
                    u_tibo_grid = u_tibo(pos_grid,:);
                    u_real_grid = u_real(pos_grid,:);
                    diff_u_grid = u_tibo_grid-u_real_grid;
                    u_err_tibo_grid = max(max(abs(diff_u_grid)));
                    H_tibo_grid = get_H(u_tibo_grid,test_type, c3);
                    H_tibo_grid_max = max(abs(H_tibo_grid));
                    if flag_plot
                        q_plot = q+2;
                        %q_plot = q;
                        p_plot = q_plot-1; %p should be smaller than q,  
                        M_plot = 2^q_plot;
                        N_plot = 2*M_plot;
                        n_plot = 2^p_plot; 
                        lambda_plot = (e-s)/n_plot;
                        x_N_plot = -b + o + lambda_plot*(0:N_plot-1);
                        x_R_plot = x_N_plot(M_plot+1:end);     
                        u_true_plot = fun_sol(x_R_plot, x_R, test_type, M, co_R, theta, deg_x, c3);
                        pos_plot = find(x_R_plot <=e & x_R_plot>=s);
                        u_tri_plot = zeros(M_plot, 3);    
                        for alpha = 1:3
                            u_R = cos_approx_engine_coef2value(coef_u(:,alpha), x_R_plot, b);
                            u_R(2:2:end) = u_R(2:2:end) -  adjustment(alpha);
                            u_tri_plot(:,alpha) = u_R;
                            u_true_this = u_true_plot(:,alpha);
                            filename_fig = ['output/ode_deg_1_system_alp_',  num2str(alpha), '_q_', num2str(q), '_t_', num2str(round(theta)), '_c_', num2str(flag_applying_constrain), '_i_', init_method ,  figtype_fig];
                            title = sprintf(char('$y_{%i, true}$ vs $y_{%i, tri}$  over [0,b]'), alpha, alpha);
                            legend_y = '$y_{true}$';
                            legend_z = '$y_{tri}$';
                            xlabel_this = '$x$';
                            ylabel_this = '$y$';
                            plot_latex_2(filename_fig, x_R_plot, u_true_this', u_R, title, legend_y, legend_z, xlabel_this, ylabel_this, plot_location )    
                        end
                    end
                    %sprintf('constrain, %i, q,  %i, deg, %i, theta, %0.1f, fval, %0.1e, u_err_tibo_grid, %0.1e,  u_rk_grid,  %0.1e, H_tri,%0.1e,  H_rk, %0.1e,\n',...
                    %flag_applying_constrain, q, deg_x, theta, fval, u_err_tibo_grid, u_rk_grid, H_tibo_grid_max, H_rk_max)
                    output(counter, :)=[flag_applying_constrain, q, deg_x, theta, fval, u_err_tibo_grid, u_rk_grid, H_tibo_grid_max, H_rk_max];   
                    counter = counter +1;
                end
            end
        end
    end
    mylib_writearray(vnames, output, ode_order_1_by_u_enhance_file);
end



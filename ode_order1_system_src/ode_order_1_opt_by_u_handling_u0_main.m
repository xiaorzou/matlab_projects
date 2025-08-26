%ode_order_1_opt_by_u_handling_u0_main.m
%{
Driver for TIBO on non-linear ODE system. 
See following model specification document:
"On Application of Trigonometric Interpolation Based Optimization solving $d$ dim Non-Linear ODE System: 
implementation specification (by $u$ with handling $u_0$)"
%}

function ode_order_1_opt_by_u_handling_u0_main()
    fourier_normalizer_const()
    global cuf_off_para 
    global flag_applying_constrains test_type weight_constrain s  e dim ps qs c3 degs init_method qqs paras ws pps
    output = zeros(length(flag_applying_constrains)*length(qqs)*length(paras)*length(degs), 9);
    counter = 1;
    ode_order_1_by_u_enhance_file = strcat('output/ode_order_1_opt_by_u_handling_u0_main_init_', init_method, '.xlsx');
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
                    co_R = co_R';
                    co_der_R = fourier_normalizer_cut_off_der1(x_R, s-delta,s,e,e+delta,cuf_off_para);
                    co_der_R = co_der_R';

                    q_2N = q+1;
                    M_2N = 2^q_2N;
                    N_2N = 2*M_2N;
                    p_2N = q_2N-1;
                    n_2N = 2^p_2N; % we should make n as large as possible,  as such, p=q-1 should be always true
                    lambda_2N = (e-s)/n_2N;
                    %m_2N = (M_2N-n_2N)/2;
                    x_2N = -b + o + lambda_2N*(0:N_2N-1);
                    x_2N_R = x_2N(M_2N+1:end);
                    co_2N_R = fourier_normalizer_cut_off(x_2N_R, s-delta,s,e,e+delta,cuf_off_para);
                    co_2N_R = co_2N_R';
                    co_2N_der_R = fourier_normalizer_cut_off_der1(x_2N_R, s-delta,s,e,e+delta,cuf_off_para);
                    co_2N_der_R = co_2N_der_R';


                    y_s_vect = fun_sol(s, x_R, test_type, M, co_R, theta, deg_x, c3); 
                    z_s_vect = fun([s, y_s_vect],x_R, test_type, M, co_R, co_der_R,theta, deg_x, ps, qs, c3);

                    u_rk = runge_kutta_ode_deg_1_variable_k(@fun_rk, x_R',y_s_vect, init_pos,x_2N_R, co_2N_R, co_2N_der_R, test_type, theta,deg_x, ps, qs, c3);
                    %z_rk = fun([x_R', u_rk]);
                    z_rk = fun([x_R', u_rk],x_R, test_type, M, co_R, co_der_R,theta, deg_x, ps, qs, c3);
                    u_rk = u_rk';
                    z_rk = z_rk';    
                    u_rk_left = zeros(M+1, dim);
                    z_rk_left = zeros(M+1, dim);
                    for i = 1:dim
                        u_rk_left(:,i) =  [0,fliplr(u_rk(i,:))]'; 
                        z_rk_left(:,i) =  [0,-fliplr(z_rk(i,:))]';
                    end
                    u_rk = u_rk';
                    z_rk = z_rk';

                    H_rk = get_H(u_rk(pos_grid,:),test_type,c3);
                    H_rk_max = max(abs(H_rk));

                    u_real = fun_sol(x_R',x_R, test_type, M, co_R, theta, deg_x, c3);
                    
                    %H_real = get_H(u_real(pos_grid,:),test_type,c3);
                    err_rk_grid = abs(u_real(pos_grid,:)-u_rk(pos_grid,:));
                    
                    err_rk_grid = max(err_rk_grid(pos_grid));
                    
                    if strcmp(test_type,'base')
                        true_fun_value = fun_sol(x_R', x_R, test_type, M, co_R, theta, deg_x, c3);
                        true_der_value = fun_der_sol(x_R', x_R, M, co_R, co_der_R, theta, deg_x, c3);
                        true_fun_value = true_fun_value';
                        true_der_value = true_der_value';
                        true_fun_value_left = zeros(M, dim);
                        true_der_value_left = zeros(M, dim);
                        for i = 1:dim
                            true_fun_value_left(:,i) =  [0,fliplr(true_fun_value(i,2:end))]'; 
                            true_der_value_left(:,i) =  [0,-fliplr(true_der_value(i,2:end))]';
                        end
                        true_fun_value = true_fun_value';

                    end
                    if strcmp(init_method, 'cheat')
                        init_guess = true_fun_value;

                    elseif strcmp(init_method, 'rk')       
                        init_guess = u_rk;
                    else
                        init_guess = zeros(M, 3);
                        for ii = 1:3
                            init_guess(:, ii) = y_s_vect(ii)*co_R;
                        end
                    end
                    r_val = get_r_value(x_R',x_R, test_type, M, co_R, co_der_R,theta, deg_x, ps, qs, c3);
                    [u_tibo, fval] = ode_order_1_opt_by_u_handling_u0_engine(@fun_w_r, @partial_u,@get_H, @get_H_der1, ...
                        y_s_vect, z_s_vect, s, e, p, init_guess, flag_applying_constrain, weight_constrain,...
                        r_val, test_type, ps, qs, c3, M, co_R);

                    u_tibo_grid = u_tibo(pos_grid,:);
                    u_real_grid = u_real(pos_grid,:);
                    diff_u_grid = u_tibo_grid-u_real_grid;
                    H_tibo_grid = get_H(u_tibo_grid,test_type, c3);
                    H_tibo_grid_max = max(abs(H_tibo_grid));
                    err_tibo = max(max(diff_u_grid));
                    output(counter, :)=[flag_applying_constrain, q, deg_x, theta, fval, err_tibo, err_rk_grid, H_tibo_grid_max, H_rk_max];   
                    counter = counter +1;                                        
                end
            end
        end           
    end
    mylib_writearray(vnames, output, ode_order_1_by_u_enhance_file); 

end

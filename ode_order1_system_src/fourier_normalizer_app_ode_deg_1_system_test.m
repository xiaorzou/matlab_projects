
%{
 test driver for general first order ODE
%}

function fourier_normalizer_app_ode_deg_1_system_test()
    addpath('D:/matlab/fourier_normalizer/fourier_normalizer_engine')  
    addpath('D:/matlab/cos_approx/cos_approx_engine')  
    addpath('D:/matlab/mylib') 
    fourier_normalizer_const()
    global cuf_off_para 

    figtype_fig = '.fig';
    init_method = 'simple';
    %init_method = 'cheat';
    %init_method = 'const';

    plot_location = 'northwest';    
    s = 1;
    e = 3*s;
    q = 7;
    p = q-1; %p should be smaller than q,  
    
    para = 2;
    theta = (0.5+para)*pi;
    
    
    M = 2^q;
    N = 2*M;
    n = 2^p; % we should make n as large as possible,  as such, p=q-1 should be always true
    lambda = (e-s)/n;
    m = (M-n)/2;
    delta = lambda*m;
    o = s - delta;
    b = e-s + 2*delta;
    
    x_N = -b + o + lambda*(0:N-1);
    x_M = x_N(1:M);
    x_R = x_N(M+1:end);
    co_R = fourier_normalizer_cut_off(x_R, s-delta,s,e,e+delta,cuf_off_para);
    co_L = [0,fliplr(co_R(2:M))];
    co_R = co_R';
  
    ps = [0.1,0.1,0.1]; 
    qs = [0.1,0.1,0.1];
    c3 = 1;

    
    dim = 3;
    
    
    % X should be always non-negative,  y1', y2', y3' are all odd function.
    function val = get_r_value(X)
        len = length(X);
        val = zeros(len, 3);
        val(:,1) = theta*cos(theta*X)- ps(1)*cos(theta*X).^2 - qs(1)*sin(theta*X);
        val(:,2) = -theta*sin(theta*X)- ps(2)*c3^2*X.^2 - qs(2)*cos(theta*X);
        val(:,3) = c3- ps(3)*sin(theta*X).^2 - c3*qs(3)*X;
    end
    
    function val = fun(XU)  % X>0 is must!  shape N , 4 X(:,1) x value,  X(:, 2) u1 value, X(:,3): u2 value, X(:,4): u3 value),  val shape: (V, 3) 
        [rows, temp] = size(XU);
        val = zeros(rows, 3);
        r_value = get_r_value(XU(:,1));
        if rows == M
            val(:,1) = (r_value(:,1) + ps(1)* XU(:,3).^2+qs(1)*XU(:,2)).*co_R;
            val(:,2) = (r_value(:,2) + ps(2)* XU(:,4).^2+qs(2)*XU(:,3)).*co_R;
            val(:,3) = (r_value(:,3) + ps(3)* XU(:,2).^2+qs(3)*XU(:,4)).*co_R;            
        else
            val(:,1) = (r_value(:,1) + ps(1)* XU(:,3).^2+qs(1)*XU(:,2));
            val(:,2) = (r_value(:,2) + ps(2)* XU(:,4).^2+qs(2)*XU(:,3));
            val(:,3) = (r_value(:,3) + ps(3)* XU(:,2).^2+qs(3)*XU(:,4));        
        end
    end
    
    function val = partial_u(XU, dim_F, dim_u)
        len = length(XU(:,1));
        u1= XU(:,2);
        u2= XU(:,3);
        u3 = XU(:,4);
        if dim_F == 1
            if dim_u == 1
                val = qs(1)*ones(len,1);
            elseif dim_u == 2
                val = 2*ps(1)*u2;
            else
                val = zeros(len,1);
            end
        elseif dim_F == 2
            if dim_u == 1
                val = zeros(len,1);
            elseif dim_u == 2
                val = qs(2)*ones(len,1);
            else
                val = 2*ps(2)*u3;
            end    
        else
            if dim_u == 1
                val = 2*ps(3)*u1;
            elseif dim_u == 2
                val = zeros(len,1);
            else
                val = qs(3)*ones(len,1);
            end  
        end
        val = val.*co_R;
    end


    function val = fun_sol(X)  %only for x>0
        len = length(X);
        val = zeros(len, 3);
        val(:,1) = sin(theta*X);
        val(:,2) = cos(theta*X);
        val(:,3) = c3*X;
    end

    % for debug purpose
    y_s_vect = fun_sol(s); 

    function val = fun_der_sol(X)  %only for x>=0
        len = length(X);
        val = zeros(len, 3);
        val(:,1) = theta*cos(theta*X);
        val(:,2) = -theta*sin(theta*X);
        val(:,3) = c3*ones(len, 1);
    end

    %u_rk = runge_kutta_ode_orde_1(@fun, x_R,y_s, m+1, []);
    %u_true = fun_sol(x_R);
    %plot(x_R, u_true, x_R, u_rk)

    true_fun_value = fun_sol(x_R);
    true_der_value = fun_der_sol(x_R);
    true_der_value_left = zeros(M, 3);
    for i = 1:3
        temp_ = true_der_value(:,i);
        true_der_value_left(:,i) = -[0; fliplr(temp_(2:M)')'];
    end
    if strcmp(init_method, 'cheat')
        init_guess = true_der_value;
        for i = 1:3
            temp_ = true_der_value(:,i).*co_R;
            init_guess(:,i) = [0; -fliplr(temp_(2:M)')'];
        end
    elseif strcmp(init_method, 'simple')
        %u_init = u_true;
        y_s_vect = fun_sol(s); 
        z_s_vect = fun([s, y_s_vect]); %
        init_guess = zeros(M,3);
        u_init = zeros(M,3);
        index_0 = m+n+1;
        init_guess(index_0,:) = -z_s_vect;   %z is odd!!!,  m+n+1 relfects to -s
        u_init(index_0,:) = y_s_vect; 
        for i = (index_0+1):M
            for ii =1:3
                u_init(i, ii) = (init_guess(i-1, ii)*lambda + u_init(i-1, ii))*co_L(i); %u is even!
            end
            %z_init(i) = fun_org([-x_M(i); u_init(i)])*co_L(i);
            init_guess(i, :) = -fun([-x_M(i), u_init(i,:)]); % z is odd,  fun is only defined over positive range!!!
            for ii = 1:3
                init_guess(i, ii) = init_guess(i, ii)*co_L(i);
            end
        end
        for i = (index_0-1):-1:2
            for ii = 1:3
                u_init(i,ii) = (init_guess(i+1, ii)*(-lambda) + u_init(i+1, ii))*co_L(i); %u is even!
            end
            init_guess(i,:) = -fun([-x_M(i), u_init(i,:)])*co_L(i); %z is odd,  fun is only defined over positive range!!!
            for ii = 1:3
                init_guess(i,ii) = init_guess(i,ii)*co_L(i); %z is odd!
            end
        end
    else
        y_s_vect = fun_sol(s); 
        z_s_vect = fun([s, y_s_vect]); %
        init_guess = zeros(M, 3);
        for ii = 1:3
            init_guess(:, ii) = -z_s_vect(ii)*co_L;
        end
    end
    
    [u_N_out, z_N_out, coef, fval] = fourier_normalizer_app_ode_deg_1_system(@fun, @partial_u, y_s_vect, s, e, p, init_guess);
    
    u_N_out = u_N_out(M+1:N,:);
    u_covered = u_N_out((m+1):(m+n+1),:);
    u_true_covered = true_fun_value((m+1):(m+n+1),:);
    u_true_der_covered = true_der_value((m+1):(m+n+1),:);
    
    z_N_out = z_N_out(M+1:N,:);
    z_covered = z_N_out((m+1):(m+n+1),:);
    
    u_err = u_covered-u_true_covered;
    z_err = z_covered - u_true_der_covered;
   
    err_grid = max(max(abs(u_err)))/max(max(abs(u_true_covered)));
    max_err_z = max(max(abs(z_err)))/max(max(abs(u_true_der_covered)));
    
    
    q_plot = 10;
    p_plot = q_plot-1; %p should be smaller than q,  
    M_plot = 2^q_plot;
    N_plot = 2*M_plot;
    n_plot = 2^p_plot; 
    lambda_plot = (e-s)/n_plot;
    %o = s - delta;
    %b = e-s + 2*delta;
    x_N_plot = -b + o + lambda_plot*(0:N_plot-1);
    
    x_R_plot = x_N_plot(M_plot+1:end);
    pos = find(x_R_plot<=e & x_R_plot>=s);  
    y_true = fun_sol(x_R_plot);
    error_plot = zeros(1,3);
    u_all = zeros(M_plot, 3);
    for alpha = 1:3
        u_R = cos_approx_engine_coef2value(coef(:,alpha), x_R_plot, b);
        u_all(:,alpha) = u_R;
        u_true = y_true(:,alpha);
        filename_fig = ['output/fourier_normalizer_app_ode_deg_1_system_test_alpha_',  num2str(alpha), '_q_', num2str(q), '_theta_', num2str(round(theta)), '_qplot_', num2str(q_plot), '_init_', init_method ,  figtype_fig];
        title = sprintf(char('$y_{%i, true}$ vs $y_{%i, opt}$  over [0,b]'), alpha, alpha);
        legend_y = '$y_{true}$';
        legend_z = '$y_{opt}$';
        xlabel_this = '$x$';
        ylabel_this = '$y$';
        plot_latex_2(filename_fig, x_R_plot, u_true', u_R, title, legend_y, legend_z, xlabel_this, ylabel_this, plot_location )
        u_R = u_R(pos); 
        u_true = u_true(pos);
        error_plot(alpha) = max(abs(u_R-u_true'))/max(abs(u_true));
    end
    error_plot = max(error_plot);
    sprintf('init_method, %s, theta, %0.1f, fval, %0.1e, err_grid, %0.1e, error_plot, %0.1e, q, %i',init_method, theta, fval, err_grid, error_plot, q)
end

%{
init_method, simple, theta, 1.6, fval, 1.8e-18,  error_plot, 4.7e-03, q, 4
init_method, simple, theta, 1.6, fval, 5.1e-18,  error_plot, 1.6e-04, q, 5
init_method, simple, theta, 1.6, fval, 1.2e-17,  error_plot, 8.9e-07, q, 6
init_method, simple, theta, 1.6, fval, 4.1e-17,  error_plot, 8.4e-09, q, 7
init_method, simple, theta, 1.6, fval, 5.0e-16,  error_plot, 2.7e-08, q, 8
init_method, const, theta, 11.0, fval, 7.5e-17,  error_plot, 1.5e-08, q, 7


init_method, simple, theta, 1.6, fval, 4.1e-17,  error_plot, 8.4e-09, q, 7
init_method, cheat, theta, 1.6, fval, 1.7e-16,   error_plot, 5.5e-08, q, 7
init_method, const, theta, 1.6, fval, 3.6e-16,   error_plot, 1.7e-08, q, 7

init_method, simple, theta, 1.6, fval, 4.1e-17,  error_plot, 8.4e-09, q, 7
init_method, simple, theta, 4.7, fval, 5.2e-17,  error_plot, 2.6e-08, q, 7
init_method, simple, theta, 7.9, fval, 7.8e-17,  error_plot, 5.6e-08, q, 7
init_method, simple, theta, 11.0, fval, 7.0e-17, error_plot, 1.6e-08, q, 7

%}



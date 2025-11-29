
%{
 test driver for general first order ODE
%}

function main_order1_singular()
    global cuf_off_para
    figtype_eps = '.eps';
    figtype_fig = '.fig';

    task = 'non_linear';
    %init_method = 'runge_kutta';
    init_method = 'simple';
    %init_method = 'cheating';
    plot_location = 'northwest';
    exclude_x0 = true;
    test_g = false;
    
    fourier_normalizer_const()
    s = 1;
    e = 3*s;
    q = 7;
    p = q-1; %p should be smaller than q,  
    
    
    M = 2^q;
    n = 2^p; % we should make n as large as possible,  as such, p=q-1 should be always true
    lambda = (e-s)/n;
    m = (M-n)/2;
    delta = lambda*m;
    o = s - delta;
    b = e-s + 2*delta;
    
    N =2*M;
    x_N = -b + o + lambda*(0:N-1);
    x_R = x_N(M+1:end);
    x_L = x_N(1:M);
    co_R = fourier_normalizer_cut_off(x_R, s-delta,s,e,e+delta,cuf_off_para);    
    co_L = [0,fliplr(co_R(2:M))];
    deg = 1;
    alpha = 0;
    theta = (0.5+alpha)*pi;
    
    function val = get_g(x) %x>=0!
        val = ones(1,length(x));
    end

    function val = get_G_R(x) %x>=0!
        if test_g
            val = get_g(x).*co_R;
        else
            val = get_g(x);
        end                
    end

    function val = get_G_L(x) %x>=0! G(x) is even 
        if test_g
            val = get_G_R(x);
            val = [0,fliplr(val(2:end))];
        else
            val = ones(1,length(x));
        end
    end    
    
    G_L = get_G_L(x_R);
    
    function val = get_r(x) %x>=0!
        val = (deg*x.^(deg-1).*cos(theta*x) - theta*x.^(deg).*sin(theta*x))-x.^(deg+1).*cos(theta*x).*(1+x.^(deg-1).*cos(theta*x));
    end
    
    function val = fun(X, r) %%x>=0! only on right side
        x = X(1,:);
        y = X(2,:);
        val = (r + x.*y+y.^2).*co_R.^2;
    end

    
    function val = fun_rk(X) %%x>=0! only on right side
        x = X(1,:);
        y = X(2,:);
        r = deg*x.^(deg-1).*cos(theta*x) - theta*x.^(deg).*sin(theta*x)-x.^(deg+1).*cos(theta*x).*(1+x.^(deg-1).*cos(theta*x));
        val = (r + x.*y+y.^2);
    end

    function val = F_for_simple(x,y, co_R_0) %%x>=0! only on right side
        r = (deg*x^(deg-1)*cos(theta*x) - theta*x^(deg)*sin(theta*x))-x^(deg+1)*cos(theta*x)*(1+x^(deg-1).*cos(theta*x));
        val = (r + x*y+y^2)*co_R_0^2;
    end


    function val = partial_u(X)
        x = X(1,:);
        y = X(2,:);
        val = (x+2*y).*co_R.^2;
    end    
    
    function val = fun_sol(x)  %only for x>=0
        val = x.^deg.*cos(theta*x);
    end
    function val = fun_der_sol(x)  %only for x>=0
        val = deg*x.^(deg-1).*cos(theta*x)-theta*x.^deg.*sin(theta*x);
    end
    
    R = get_r(x_R);   
    y_s = fun_sol(s);
    z_s = fun_der_sol(s);
    u_real = fun_sol(x_R);
    u_b = fun_sol(b);
    z_real = fun_der_sol(x_R);
    u_real_L = [u_b,fliplr(u_real(2:M))];
    z_real_L = [0,fliplr(z_real(2:M))];    
    
    z_M = [0,-fliplr(z_real(2:M))];
    grid_pos = [m+1:m+n+1]; 
    
    eta_test = zeros(1,M);
    eta_test(grid_pos) = 1; %we should refine this so that it have more weights on singular point. 
    eta_test = [0,fliplr(eta_test(2:M))];
    if test_g 
        %eta = eta_test;
        eta = ones(1,M);
    else
        eta = ones(1,M);
    end
        
    
    %verify objective function!
    xu_R = [x_R;u_real];
    F_R = fun(xu_R, R); 
    F_M = [0,-fliplr(F_R(2:M))];
    myobj = max(abs(eta_test.*(G_L.*z_M-F_M)));
    sprintf('myobj is %0.1e', myobj)

    if strcmp(init_method, 'simple')
       %u_init = u_true;
        z_init = zeros(1,M);
        u_init = zeros(1,M);
        index_0 = m+n+1;
        z_init(index_0) = -z_s;  %z is odd, we are working at x=-s
        u_init(index_0) = y_s;   %u is even , we are working at x=-s
        for i = (index_0+1):M
            u_init(i) = (z_init(i-1)*lambda + u_init(i-1));
             if abs(G_L(i)) ~= 0
                %z_init(i) = -fun_rk([-x_L(i); u_init(i)])*co_L(i).^2/G_L(i); %z is odd!
                z_init(i) = -F_for_simple(-x_L(i),u_init(i), co_L(i))/G_L(i);
             else
                 z_init(i) = 0;
             end
        end
        for i = (index_0-1):-1:1
            u_init(i) = (z_init(i+1)*(-lambda) + u_init(i+1)); %u is even!
            if abs(G_L(i))~=0
                %z_init(i) = -fun_rk([-x_L(i); u_init(i)])*co_L(i).^2/G_L(i); %z is odd!
                z_init(i) = -F_for_simple(-x_L(i),u_init(i), co_L(i))/G_L(i); %z is odd!
            else
                z_init(i) = 0;
            end
        end
    elseif strcmp(init_method, 'cheating')
        z_init = [0,-fliplr(z_real(2:M))];       
    elseif strcmp(init_method, 'rk4')
        u_rk = runge_kutta_ode_orde_1(fun_rk, x_R,y_s, m+1, []);
        z_rk =  fun([x_R; u_rk]);
        z_init =  [0,-fliplr(z_rk(2:M))];
    end
    if exclude_x0
        z_init =  [z_init(2:m+n),z_init(m+n+2:end)];
    else
        z_init =  [z_init(1:m+n),z_init(m+n+2:end)];
    end
    [z_N, u_N, opt_error] = engine_order1_singular(@fun, @partial_u, y_s, z_s, s, e, p,q, G_L, R,eta, z_init, exclude_x0);
    pos = find(x_N >= s & x_N<=e);
    x_M_plot = x_N(pos);
    u_M_plot = u_N(pos);

    u_true = fun_sol(x_M_plot);
    u_der = fun_der_sol(x_M_plot);
    

    

    %handle error
    error_fig_file = ['output/sect7_fig5_', task, '_q_', num2str(q), '_deg_', num2str(deg), '_theta_', num2str(round(theta)), '_init_', init_method ,  figtype_fig];
    title = strcat('The comparison of the errors between $intp$ and $benc$');
    legend_y = '$intp$';
    legend_z = '$benc$';
    xlabel_this = 'x';
    ylabel_this = 'y';
    u_rk = runge_kutta_ode_orde_1(@fun_rk, x_M_plot,y_s, 1, []);
    u_rk_invest = runge_kutta_ode_orde_1(@fun_rk, x_M_plot,y_s, 1, u_true);
    error_opt = u_M_plot - u_true;
    error_rk = u_rk - u_true;
    error_rk_invest = u_rk_invest - u_true;
    plot_latex_2(error_fig_file, x_N(pos), error_opt,  error_rk_invest,title, legend_y, legend_z, xlabel_this, ylabel_this,plot_location) 
    
    
    max_error_q_opt = [max(abs(error_opt)), max(abs(error_rk)), max(abs(error_rk_invest))];
    
    %handle IA
    max_der = max(abs(u_der));
    max_error_rk = max(abs(error_rk));
    u_der = max_error_rk*u_der/max_der;
    title = strcat('The derivative impact on rk4 error');
    legend_y = '$deri$';
    legend_z = '$rk4$';
    ia_der_eps_file = ['output/fourier_normalizer_app_ode_deg_1_ia_der_init_', task, '_q_', num2str(q), '_deg_', num2str(deg), '_theta_', num2str(round(theta)), '_init_', init_method ,  figtype_eps];
    ia_der_fig_file = ['output/fourier_normalizer_app_ode_deg_1_ia_der_init_', task, '_q_', num2str(q), '_deg_', num2str(deg), '_theta_', num2str(round(theta)), '_init_', init_method ,  figtype_fig];
    plot_latex_2(ia_der_eps_file, x_N(pos), u_der,  error_rk, title, legend_y, legend_z, xlabel_this, ylabel_this,plot_location) 
    plot_latex_2(ia_der_fig_file, x_N(pos), u_der,  error_rk, title, legend_y, legend_z, xlabel_this, ylabel_this,plot_location) 
    

    %handle difference
    x_shift = x_N(pos);
    x_shift = x_shift(2:end);
    diff_true = u_true(2:end)-u_true(1:end-1);
    diff_opt = u_M_plot(2:end)-u_M_plot(1:end-1)-diff_true;
    diff_rk = u_rk(2:end)-u_rk(1:end-1)-diff_true;
    diff_rk_invest = u_rk_invest(2:end)-u_rk_invest(1:end-1)-diff_true;
    title = strcat('The comparison of the differences among $intp$, $rk4$ and $benc$');
    legend_y = '$intp$';
    legend_z = '$rk4$';
    legend_w = '$benc$';
    ia_der_eps_file = ['output/sect7_fig6_', task, '_q_', num2str(q), '_deg_', num2str(deg), '_theta_', num2str(round(theta)), '_init_', init_method ,  figtype_eps];
    ia_der_fig_file = ['output/sect7_fig6_', task, '_q_', num2str(q), '_deg_', num2str(deg), '_theta_', num2str(round(theta)), '_init_', init_method ,  figtype_fig];
    plot_latex_3(ia_der_eps_file, x_shift, diff_opt,  diff_rk, diff_rk_invest, title, legend_y, legend_z, legend_w, xlabel_this, ylabel_this,plot_location) 
    plot_latex_3(ia_der_fig_file, x_shift, diff_opt,  diff_rk, diff_rk_invest, title, legend_y, legend_z, legend_w, xlabel_this, ylabel_this,plot_location) 

    
    

    %handle y and u
    [coef, ~] = cos_approx_engine_value2coef_wo_loop(u_N);
    q_plot = 10;
    p = q_plot-1; %p should be smaller than q,  
    M = 2^q_plot;
    N = 2*M;
    n = 2^p; 
    lambda = (e-s)/n;
    %o = s - delta;
    %b = e-s + 2*delta;
    x_N = -b + o + lambda*(0:N-1);
    
    x_R = x_N(M+1:end);
    x_R_selected = find(x_R<=e & x_R>=s);    
    u_R = cos_approx_engine_coef2value(coef, x_R, b);
  
    y_true = fun_sol(x_R);
    max_error_q_plot = max(abs(y_true(x_R_selected) - u_R(x_R_selected)));

    filename_eps = ['output/sect7_fig7_', task, '_q_', num2str(q), '_deg_', num2str(deg), '_theta_', num2str(round(theta)), '_qplot_', num2str(q_plot), '_init_', init_method ,  figtype_eps];
    filename_fig = ['output/sect7_fig7_', task, '_q_', num2str(q), '_deg_', num2str(deg), '_theta_', num2str(round(theta)), '_qplot_', num2str(q_plot), '_init_', init_method ,  figtype_fig];
    title = strcat('The plots of $y$ and $u$ over [0,b]');
    legend_y = '$y$';
    legend_z = '$u$';
    xlabel_this = '$x$';
    ylabel_this = '$y$ and $u$';
    plot_latex_2(filename_eps, x_R, y_true, u_R, title, legend_y, legend_z, xlabel_this, ylabel_this, plot_location)
    plot_latex_2(filename_fig, x_R, y_true, u_R, title, legend_y, legend_z, xlabel_this, ylabel_this, plot_location )

    
    %handle performance results and parameters
    fourier_normalizer_app_ode_deg_1_general_max_error_file = ['output/sect7_tab6_', task, '_q_', num2str(q), '_deg_', num2str(deg), '_theta_', num2str(round(theta)), '_init_', init_method ,  '.xlsx'];
    vnames = {'intp', 'rk4', 'benc', 'plot', 'opt', 'q'};
    
    res = zeros(1,6);
    res(1,1)= max_error_q_opt(1);
    res(1,2)= max_error_q_opt(2);
    res(1,3)= max_error_q_opt(3);
    res(1,4)= max_error_q_plot;
    res(1,5)= opt_error;
    res(1,6)= q;
    mylib_writearray(vnames, res, fourier_normalizer_app_ode_deg_1_general_max_error_file);
    
end

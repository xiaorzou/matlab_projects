
%{
 test driver for general first order ODE
%}

function main_ode_order1()
    figtype_eps = '.eps';
    figtype_fig = '.fig';

    task = 'non_linear';
    %init_method = 'runge_kutta';
    init_method = 'simple';
    plot_location = 'northwest';
    
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
    
    if strcmp(task, 'linear')
        y_s = 2;
    elseif strcmp(task, 'non_linear')
        y_s = 0;
    end
    
    deg = 1;
    alpha = 1;
    theta = (0.5+alpha)*pi;
    
    function val = fun(X)
        if strcmp(task, 'linear')
            val = (abs(X(1,:)).^2 -abs(X(1,:)).^2.*X(2,:));
        elseif strcmp(task, 'non_linear')
            x = abs(X(1,:));
            y = X(2,:);
            f = deg*x.^(deg-1).*cos(theta*x) - theta*x.^(deg).*sin(theta*x)-x.^(deg+1).*cos(theta*x).*(1+x.^(deg-1).*cos(theta*x));
            val = f + x.*y+y.^2;
        else
            disp(['not implemented with task ', task])
            return;
        end
    end

    z_s = fun([s; y_s]); %old extension,  it is valued at x=-s
    function val = partial_u(X)
        if strcmp(task, 'linear')
            val = -abs(X(1,:)).^2;
        elseif strcmp(task, 'non_linear')
            x = abs(X(1,:));
            y = X(2,:);
            val = x+2*y;
        else
            disp(['not implemented with task ', task])
            return;
        end
    end

    function val = fun_sol(X)  %only for x<=0
        if strcmp(task, 'linear')
            x = abs(X);
            val = ((y_s-1)*exp(-(x.^3-s^3)/3) +1); %only for x<=0
        elseif strcmp(task, 'non_linear')
            x = abs(X);
            val = x.^deg.*cos(theta*x);
        else
            disp(['not implemented with task ', task])
            return;
        end
    end

    function val = fun_der_sol(X)  %only for x<=0
        if strcmp(task, 'linear')
            x = abs(X);
            val = (y_s-1)*exp(-(x.^3-s^3)/3).*(-x.^2);
        elseif strcmp(task, 'non_linear')
            x = abs(X);
            val = deg*x.^(deg-1).*cos(theta*x)-theta*x.^deg.*sin(theta*x);
        else
            disp(['not implemented with task ', task])
            return;
        end
    end
    
    [x_N, u_N, u_init_N, z_N, z_init_N, pos, opt_error] = engine_ode_order1(@fun, @partial_u, y_s, z_s, s, e, p,q, init_method);
    x_M_plot = x_N(pos);
    u_M_plot = u_N(pos);
    x_M_plot_size = length(x_M_plot);
    u_true = fun_sol(x_M_plot);
    u_der = fun_der_sol(x_M_plot);

    %handle error
    error_eps_file = ['output/sect7_fig5_', task, '_q_', num2str(q), '_deg_', num2str(deg), '_theta_', num2str(round(theta)), '_init_', init_method ,  figtype_eps];
    error_fig_file = ['output/sect7_fig5_', task, '_q_', num2str(q), '_deg_', num2str(deg), '_theta_', num2str(round(theta)), '_init_', init_method ,  figtype_fig];
    title = strcat('The comparison of the errors between $intp$ and $benc$');
    legend_y = '$intp$';
    legend_z = '$benc$';
    xlabel_this = 'x';
    ylabel_this = 'y';
    u_rk = runge_kutta_ode_orde_1(@fun, x_M_plot,y_s, 1, []);
    u_rk_invest = runge_kutta_ode_orde_1(@fun, x_M_plot,y_s, 1, u_true);
    error_opt = u_M_plot - u_true;
    error_rk = u_rk - u_true;
    error_rk_invest = u_rk_invest - u_true;
    plot_latex_2(error_eps_file, x_N(pos), error_opt,  error_rk_invest,title, legend_y, legend_z, xlabel_this, ylabel_this,plot_location) 
    plot_latex_2(error_fig_file, x_N(pos), error_opt,  error_rk_invest,title, legend_y, legend_z, xlabel_this, ylabel_this,plot_location) 
    
    %{
    error_output = zeros(x_M_plot_size, 3);
    error_output(:,1) = x_M_plot;
    error_output(:,2) = error_opt;
    %error_output(:,3) = error_rk;
    error_output(:,3) = error_rk_invest;
    vnames = {'x','intp','benc'};   
    fourier_normalizer_app_ode_deg_1_general_error_xlsx_file = ['output/fourier_normalizer_app_ode_deg_1_error_', task, '_q_', num2str(q) , '_deg_', num2str(deg), '_theta_', num2str(round(theta)), '_init_', init_method , '.xlsx'];
    mylib_writearray(vnames, error_output, fourier_normalizer_app_ode_deg_1_general_error_xlsx_file);
    %}
    
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
    
    %fourier_normalizer_app_ode_deg_1_general_ia_der_xlsx_file = ['output/fourier_normalizer_app_ode_deg_1_ia_der_', task, '_q_', num2str(q) , '_deg_', num2str(deg), '_theta_', num2str(round(theta)), '_init_', init_method , '.xlsx'];
    %{
    ia_der_output = zeros(x_M_plot_size, 3);
    ia_der_output(:,1) = x_N(pos);
    ia_der_output(:,2) = u_der;
    ia_der_output(:,3) = error_rk;
    vnames = {'x','deri','error'};   
    mylib_writearray(vnames, ia_der_output, fourier_normalizer_app_ode_deg_1_general_ia_der_xlsx_file);
    %}
    

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
    %{
    diff_output = zeros(length(x_shift), 4);
    diff_output(:,1) = x_shift;
    diff_output(:,2) = diff_opt;
    diff_output(:,3) = diff_rk;
    diff_output(:,4) = diff_rk_invest;  
    vnames = {'x','intp','rk4', 'benc'};   
    mylib_writearray(vnames, diff_output, fourier_normalizer_app_ode_deg_1_general_ia_der_xlsx_file);
   %}
    
    

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
    %{
    fourier_normalizer_app_ode_deg_1_general_yandu_xlsx_file = ['output/fourier_normalizer_app_ode_deg_1_yandu_', task, '_q_', num2str(q) , '_deg_', num2str(deg), '_theta_', num2str(round(theta)), '_init_', init_method , '.xlsx'];
    yandu_output = zeros(length(x_R), 3);
    yandu_output(:,1) = x_R;
    yandu_output(:,2) = y_true;
    yandu_output(:,3) = u_R;
    vnames = {'x','y','u'};   
    mylib_writearray(vnames, yandu_output, fourier_normalizer_app_ode_deg_1_general_yandu_xlsx_file);
    %}
    
    %handle performance results and parameters
    fourier_normalizer_app_ode_deg_1_general_max_error_file = ['output/sect7_tab6_', task, '_q_', num2str(q), '_deg_', num2str(deg), '_theta_', num2str(round(theta)), '_init_', init_method ,  '.xlsx'];
    vnames = {'intp', 'rk4', 'benc', 'plot', 'opt', 'lambda', 'max_der'};
    
    res = zeros(1,7);
    res(1,1)= max_error_q_opt(1);
    res(1,2)= max_error_q_opt(2);
    res(1,3)= max_error_q_opt(3);
    res(1,4)= max_error_q_plot;
    res(1,5)= opt_error;
    res(1,6)= lambda;
    res(1,7)= max_der;
    mylib_writearray(vnames, res, fourier_normalizer_app_ode_deg_1_general_max_error_file);
    
end

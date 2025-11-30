
% test driver for results of 4.2 of the following paper
% 
% "Trigonometric Interpolation on Non-Periodic Functions and its Applications"
% 
% to get the results for theta = pi/2, set alpha =0,  
% to get the results for theta = 3*pi/2, set alpha =1  
% 

function main_ode_order1()
    task = 'non_linear';
    figtype_fig = '.fig';


    %init_method = 'runge_kutta';
    init_method = 'simple';
    plot_location = 'northwest';
    
    smoothing = true;
    save_fig = true;
    
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
    alpha = 0;
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
    
    [x_N, u_N, u_init_N, z_N, z_init_N, pos, opt_error] = engine_ode_order1(@fun, @partial_u, y_s, z_s, s, e, p,q, init_method, smoothing);
    x_M_plot = x_N(pos);
    u_M_plot = u_N(pos);
    u_true = fun_sol(x_M_plot);
    u_der = fun_der_sol(x_M_plot);

    %handle error

    u_rk = runge_kutta_ode_orde_1(@fun, x_M_plot,y_s, 1, []);
    u_rk_invest = runge_kutta_ode_orde_1(@fun, x_M_plot,y_s, 1, u_true);
    error_opt = u_M_plot - u_true;
    error_rk = u_rk - u_true;
    error_rk_invest = u_rk_invest - u_true;
    if save_fig        
        error_fig_file = ['output/fig5_', task, '_q_', num2str(q), '_deg_', num2str(deg), '_theta_', num2str(round(theta)), '_init_', init_method ,  figtype_fig];
        title = strcat('The comparison of the errors between $intp$ and $benc$');
        legend_y = '$intp$';
        legend_z = '$benc$';
        xlabel_this = 'x';
        ylabel_this = 'y';
        plot_latex_2(error_fig_file, x_N(pos), error_opt,  error_rk_invest,title, legend_y, legend_z, xlabel_this, ylabel_this,plot_location) 
    end
    
    max_error_q_opt = [max(abs(error_opt)), max(abs(error_rk)), max(abs(error_rk_invest))];
    %handle IA
    max_der = max(abs(u_der));
    
    if save_fig        
        max_error_rk = max(abs(error_rk));
        u_der = max_error_rk*u_der/max_der;
        title = strcat('The derivative impact on rk4 error');
        legend_y = '$deri$';
        legend_z = '$rk4$';
        ia_der_fig_file = ['output/fourier_normalizer_app_ode_deg_1_ia_der_init_', task, '_q_', num2str(q), '_deg_', num2str(deg), '_theta_', num2str(round(theta)), '_init_', init_method ,  figtype_fig];
        plot_latex_2(ia_der_fig_file, x_N(pos), u_der,  error_rk, title, legend_y, legend_z, xlabel_this, ylabel_this,plot_location) 
    end
    
    %handle difference

    if save_fig
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
   
        ia_der_fig_file = ['output/fig5_', task, '_q_', num2str(q), '_deg_', num2str(deg), '_theta_', num2str(round(theta)), '_init_', init_method ,  figtype_fig];
        plot_latex_3(ia_der_fig_file, x_shift, diff_opt,  diff_rk, diff_rk_invest, title, legend_y, legend_z, legend_w, xlabel_this, ylabel_this,plot_location) 
    end
    

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
    if save_fig
        
        filename_fig = ['output/fig6_', task, '_q_', num2str(q), '_deg_', num2str(deg), '_theta_', num2str(round(theta)), '_qplot_', num2str(q_plot), '_init_', init_method ,  figtype_fig];
        title = strcat('The plots of $y$ and $u$ over [0,b]');
        legend_y = '$y$';
        legend_z = '$u$';
        xlabel_this = '$x$';
        ylabel_this = '$y$ and $u$';
        
        plot_latex_2(filename_fig, x_R, y_true, u_R, title, legend_y, legend_z, xlabel_this, ylabel_this, plot_location )
    end
    %handle performance results and parameters
    
    if smoothing
        fourier_normalizer_app_ode_deg_1_general_max_error_file = ['output/tab6_', task, '_q_', num2str(q), '_deg_', num2str(deg), '_theta_', num2str(round(theta)), '_init_', init_method ,  '.xlsx'];
    else
        fourier_normalizer_app_ode_deg_1_general_max_error_file = ['output/tab6_', task, '_q_', num2str(q), '_deg_', num2str(deg), '_theta_', num2str(round(theta)), '_init_', init_method ,  '_no_smoothing.xlsx'];
    end        
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

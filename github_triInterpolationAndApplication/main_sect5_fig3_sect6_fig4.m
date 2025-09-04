%
% the function generate Figure 3 in section 5, Figure 4 in section 6 of the following paper
% "On Trigonometric Interpolation and Its Applications II"
% 
%

function  main_sect5_fig3_sect6_fig4()  
    location = 'northwest';     
    % default setting used in the paper
    s = 2;
    e = 3; 
    q = 7;
    %end of default setting
    M = 2^q;
    p = q-1;
    n = 2^p;
    
    lambda = (e-s)/n;
    m = (M-n)/2;
    delta = lambda*m;
    l = s - delta;
    
    %model  = 'default'; % used in paper
    %model = 'exp'; % used in paper
    %model = 'cos_alpha'; % used in paper
    %model = 'sin_alpha'; % used in paper
    %model = 'ln'; % not used in paper
    %model  = 'power'; % not used in paper
    
    models = {'default', 'exp', 'cos_alpha', 'sin_alpha'};
    for model = models
        if strcmp(model, 'default')
            fun_name = '$(x-2.5)^2$';
            color_code = 'b';
        elseif strcmp(model, 'power')
            fun_name = '$(1-(x/2)^2)^2$';
            color_code = 'r';
        elseif strcmp(model, 'cos_alpha')
            fun_name = '$\cos(\pi x)$';
            color_code = 'g';
        elseif strcmp(model, 'sin_alpha')
            fun_name = '$\sin(\pi x)$';
            color_code = 'c';
        elseif strcmp(model, 'exp')    
            fun_name = '$\exp(-x^2)$';
            color_code = 'm';
        elseif strcmp(model, 'ln')    
            fun_name = '$\log^2(x^2+1)/log^2(5)$';
            color_code = 'k';
        end
        
        if strcmp(model, 'default')
            fun = @(x) (x-2.5).^2;
            para = 2;
        else 
            para = 2;
            o = s-delta;
            b = e+delta-o;
            fun = functor_fun(model, para, b);
        end

        coef_co = fourier_normalizer_get_coef_by_fun(fun, s,e,p,q);

        %{
        coef_n_impact = coef_co';
        file_name_coef_n_impact = ['output/test_on_consecutive_error_s_',  num2str(s), '_e_', num2str(e), '_p_', num2str(p), '_q_', num2str(q), '_model_', model, '.xlsx'];
        vnames = {'n1'};
        mylib_writearray(vnames, coef_n_impact, file_name_coef_n_impact);
        %}
        plot_size = 2^12;
        x_plot = l + (2*delta+(e-s))/plot_size*(0:plot_size);
        y_plot = fourier_normalizer_coef2value_by_sepq(coef_co, s, e, p,q, x_plot);
        y_true = fun(x_plot);
        z = find (x_plot>=2 & x_plot <=3);
        fig_type = '.fig';
        legend_y = ['$f$'];
        legend_z = ['$\hat{f}_M$'];
        xlabel_this = 'x';
        ylabel_this = 'y';
        %{
        if strcmp(model, 'default')
            title = strcat('$f(x)=(x-2.5)^2$ vs $\hat{f}_M(x)$ over $[s,e]$');
            filename_e2s = ['output/test_on_consecutive_error_q_', num2str(q), fig_type];
        else
            title = strcat('f(x)=', fun_name,' vs $\hat{f}_M(x)$ over $[s,e]$');
            filename_e2s = ['output/hf_sepq_model_', model, num2str(round(para)) , num2str(q),  fig_type];
        end
        %plot_latex_2(filename_e2s, x_plot(z), y_true(z), y_plot(z), title, legend_y, legend_z, xlabel_this, ylabel_this,location ); 
        %}
        max_this = max(abs(y_true(z)));

        plot_skip = 2^0;
        y_true_this = y_true(z)/max_this;
        y_plot_this = y_plot(z)/max_this;
        y_true_this = y_true_this(1:plot_skip:end);
        y_plot_this = y_plot_this(1:plot_skip:end);
        x_plot_this = x_plot(z);
        x_plot_this = x_plot_this(1:plot_skip:end);

        error = y_plot_this - y_true_this;

        x_plot_this = x_plot_this(2:end);
        diff_error = error(2:end)-error(1:end-1);
        close all
        title_this = strcat('$diff_{err}$ with f=', fun_name);
        legend_this = '$diff_{err}$';
        file_name = ['output/Fig_4_sect6_', model{1}, num2str(round(para)), num2str(q), fig_type];
        plot_latex(file_name, x_plot_this, diff_error, title_this, legend_this, xlabel_this, ylabel_this, location, color_code);

        if strcmp(model, 'default')
            filename = ['output/Fig_3_sect5_', num2str(q), fig_type];
            title = strcat('$f(x)=(x-2.5)^2$ vs $\hat{f}_M(x)$ over $[s-\delta,e+\delta]$');
            plot_latex_2(filename, x_plot, y_true, y_plot, title, legend_y, legend_z, xlabel_this, ylabel_this,location ); 
        end
    end
end

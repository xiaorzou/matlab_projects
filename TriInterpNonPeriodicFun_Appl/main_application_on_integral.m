%
% the function generate table 3  and Table 4 of the following paper
% "Trigonometric Interpolation on Non-Periodic Functions and its Applications"
%

function main_application_on_integral()
    m_powers = [8];
    plot_size = 2^12;
    flag_types = {'power', 'cos'};
    fig_type = '.fig';
    flag_plot = false;
    s = -1;
    e = 1; 
    %deltas = [1*e, 2*e];
    deltas = [1*e];
    %cut_off_paras = [0.5,1,5,10];
    cut_off_paras = [0.5];    
    for flag_type = flag_types
        if strcmp(flag_type, 'power')
            paras = [4,8,10];
        elseif strcmp(flag_type, 'cos')
            paras = [ 1,10,100];
        end
        counter = 1;
        output = zeros(length(m_powers)*length(paras), 5+4+1);
        x_plot =  fourier_normalier_get_grid(plot_size, 0, e);
        output_int = zeros(length(m_powers)*length(paras), 2+2+2+2);
        for cut_off_para = cut_off_paras
            for delta = deltas
                right = e + delta;
                left = -right;
                for m_power = m_powers
                    q = m_power;
                    term = 2^m_power;
                    M = term;
                    N = 2*M;
                    x_e = (0:(N-1))*2*e/N-e; 
                    x_right = (0:(N-1))*2*right/N-right; 
                    for para = paras
                        f_e = getf(x_e, e, para, flag_type);
                        %f_a = (1-(x_a/a).^2).^para;
                        [coef, ~] = cos_approx_engine_value2coef_wo_loop(f_e);
                        f_cos_plot = cos_approx_engine_coef2value(coef, x_plot, e);
                        f_cos_plot_der1 = cos_approx_engine_coef2value_enhance(coef, x_plot, e, 1); 
                        f_cos_plot_der2 = cos_approx_engine_coef2value_enhance(coef, x_plot, e, 2); 

                        %h = cos_approx_engine_coef2value(h_coef_e, x_b, b);
                        h = fourier_normalizer_cut_off(x_right, left, s,e,right,cut_off_para);
                        f_right = getf_true(x_right, para,flag_type); 
                        fh = h.*f_right;
                        %fh = h.*(1-(x_b/a).^2).^para;
                        [coef_fh, ~] = cos_approx_engine_value2coef_wo_loop(fh);

                        x_trapz = (0:(N))*2*e/N-e; 
                        y = getf_true(x_trapz, para, flag_type);

                        true_integration = get_true_integration(e,para,flag_type);
                        output_int(counter, :) = [cut_off_para, right, q, para,  log10(abs(getInt(coef_fh, e,right, M)-true_integration)), log10(abs(trapz(x_trapz,y)-true_integration)), log10(abs(mysimpson(x_trapz,y)-true_integration)), true_integration];
                        fh_cos_plot = cos_approx_engine_coef2value(coef_fh, x_plot, right);

                        fh_cos_plot_der1 = cos_approx_engine_coef2value_enhance(coef_fh, x_plot, right, 1);
                        fh_cos_plot_der2 = cos_approx_engine_coef2value_enhance(coef_fh, x_plot, right, 2);

                        f_plot = getf(x_plot, e, para, flag_type);
                        f_plot_der1 = getf_der1(x_plot, e, para, flag_type);
                        f_plot_der2 = getf_der2(x_plot, e, para, flag_type);


                        f_plot_max = max(abs(f_plot));
                        f_plot_max_der1 = max(abs(f_plot_der1));
                        f_plot_max_der2 = max(abs(f_plot_der2));
                        output(counter,:) = [cut_off_para, right, m_power, para, log10(max(abs(f_plot-f_cos_plot))/f_plot_max), log10(max(abs(f_plot-fh_cos_plot))/f_plot_max), ...
                            log10(max(abs(f_plot_der1-f_cos_plot_der1))/f_plot_max_der1), log10(max(abs(f_plot_der1-fh_cos_plot_der1))/f_plot_max_der1), ...
                            log(max(abs(f_plot_der2-f_cos_plot_der2))/f_plot_max_der2), log10(max(abs(f_plot_der2-fh_cos_plot_der2))/f_plot_max_der1)];
                        if flag_plot && delta == deltas(1)
                            filename = ['output/test_on_integration_', flag_type, '_q_', num2str(q), '_para_', num2str(para), '_cut_off_' ,fig_type];
                            title = strcat('$y(x)=x^n, \quad n=', num2str(para), ',q=', num2str(q), '$');
                            legend_y = ['$y(x)$'];
                            legend_z = ['$y_{cos}(x)$'];
                            legend_w = ['$yh_{cos}(x)$'];
                            xlabel_this = 'x';
                            ylabel_this = '$y$';
                            s_p = 1;
                            e_p = 75;
                            plot_latex_3(filename, x_plot(s_p:e_p), f_plot(s_p:e_p), f_cos_plot(s_p:e_p),fh_cos_plot(s_p:e_p), title, legend_y, legend_z,legend_w, xlabel_this, ylabel_this )
                            filename_der1_yz = ['output/test_on_integration_der1_yz_', flag_type, '_q_', num2str(q), '_para_', num2str(para), '_cut_off_',fig_type];
                            filename_der1_yw = ['output/test_on_integration_der1_yw_', flag_type, '_q_', num2str(q), '_para_', num2str(para), '_cut_off_',fig_type];
                            title = strcat('$y^{(1)}(x)=n x^{n-1}, \quad n=', num2str(para), ',q=', num2str(q), '$');
                            legend_y = ['$y^{(1)}(x)$'];
                            legend_z = ['$y^{(1)}_{cos}(x)$'];
                            legend_w = ['$yh^{(1)}_{cos}(x)$'];
                            xlabel_this = 'x';
                            ylabel_this = '$y^{(1)}$';
                            s_p = 1;
                            plot_latex_2(filename_der1_yz, x_plot(s_p:e_p), f_plot_der1(s_p:e_p), f_cos_plot_der1(s_p:e_p), title, legend_y, legend_z, xlabel_this, ylabel_this )  
                            plot_latex_2(filename_der1_yw, x_plot(s_p:e_p), f_plot_der1(s_p:e_p), fh_cos_plot_der1(s_p:e_p), title, legend_y, legend_w, xlabel_this, ylabel_this )  

                            filename_der2_yz = ['output/test_on_integration_der2_yz_', flag_type, '_q_', num2str(q), '_para_' num2str(para), '_cut_off_', num2str(n) ,fig_type];
                            filename_der2_yw = ['output/test_on_integration_der2_yw_', flag_type, '_q_', num2str(q), '_para_' num2str(para), '_cut_off_', num2str(n) ,fig_type];
                            title = strcat('$y^{(2)}(x)=n(n-1)x^{n-2}, \quad n=', num2str(para), ',q=', num2str(q), '$');
                            legend_y = ['$y^{(2)}y(x)$'];
                            legend_z = ['$y^{(2)}_{cos}(x)$'];
                            legend_w = ['$yh^{(2)}_{cos}(x)$'];
                            xlabel_this = 'x';
                            ylabel_this = '$y^{(2)}$';
                            s_p = 1;
                            plot_latex_2(filename_der2_yz, x_plot(s_p:e_p), f_plot_der2(s_p:e_p), f_cos_plot_der2(s_p:e_p), title, legend_y, legend_z, xlabel_this, ylabel_this ) 
                            plot_latex_2(filename_der2_yw, x_plot(s_p:e_p), f_plot_der2(s_p:e_p), fh_cos_plot_der2(s_p:e_p), title, legend_y, legend_w, xlabel_this, ylabel_this )
                        end             
                        counter = counter+1;
                    end
                end
            end
        end
        fourier_normalizer_test_file = ['output/table_3_', flag_type{1}, '.xlsx'];
        vnames = {'cut_off_para','right', 'q', 'para','max_error_f', 'max_error_hf','max_error_f_der1', 'max_error_hf_der1','max_error_f_der2', 'max_error_hf_der2'};
        mylib_writearray(vnames, output, fourier_normalizer_test_file)
        fourier_normalizer_test_int_file = ['output/table_4_5_', flag_type{1}, '.xlsx'];
        vnames = {'cut_off_para','right','q', 'para', 'int_cos', 'int_trapz', 'int_simplson', 'int_true'};
        mylib_writearray(vnames, output_int, fourier_normalizer_test_int_file)
    end
end


function val = get_true_integration(e,para,flag_type)
    if strcmp(flag_type, 'cos')
        val = 2/(para)*sin(para*e);
    elseif strcmp(flag_type, 'power')
        val = 2*e/(para+1);
    end
end



function val = getf_true_cos(x, para)
    val = cos(para*x);
end

function val = getf_true_power(x, deg)
    val = x.^(deg);
end

function val = getf_true(x, para, flag_type)
    if strcmp(flag_type, 'cos')
        val = getf_true_cos(x, para);
    elseif strcmp(flag_type, 'power')
        val = getf_true_power(x, para);
    end
end


function val = getf_cos(x, a, para)
    N = length(x);
    val = zeros(1,N);
    for i =1:N
        if x(i)<=a && x(i)>=-a
            val(i) = cos(para*x(i));
        elseif x(i)<-a
            val(i) = cos(para*(-2*a-x(i))); %(1-((-2*a-x(i))/a)^2)^para;
        else
            val(i) = cos(para*(2*a-x(i)));%(1-((2*a-x(i))/a)^2)^para;
        end
    end
end


function val = getf_power(x, e, deg)
    N = length(x);
    val = zeros(1,N);
    for i =1:N
        if x(i)<=e && x(i)>=-e
            val(i) = x(i)^(deg);
        elseif x(i)<-e
            val(i) = (-2*e-x(i))^(deg);
        else
            val(i) = (2*e-x(i))^(deg);
        end
    end
end

function val = getf(x, e, para, flag_type)
    if strcmp(flag_type, 'cos')
        val = getf_cos(x, e, para);
    elseif strcmp(flag_type, 'power')
        val = getf_power(x, e, para);
    end
end

function val = getf_der1_cos(x, a, para)
    N = length(x);
    val = zeros(1,N);
    for i =1:N
        if x(i)<=a && x(i)>=-a
            val(i) = -para* sin(para*x(i)); %para*(1-(x(i)/a)^2)^(para-1)*(-2)*x(i)/a^2;
        elseif x(i)<-a
            val(i) = para* sin(para*(-2*a-x(i))); %para*(1-((-2*a-x(i))/a)^2)^(para-1)*(2)*(-2*a-x(i))/a^2;
        else
            val(i) = para* sin(para*(2*a-x(i))); %para*(1-((2*a-x(i))/a)^2)^(para-1)*(2)*(2*a-x(i))/a^2;
        end
    end
end

function val = getf_der1_power(x, e, deg)
    N = length(x);
    val = zeros(1,N);
    for i =1:N
        if x(i)<=e && x(i)>=-e
            val(i) = deg*x(i)^(deg-1);
        elseif x(i)<-e
            val(i) = -deg*(-2*e-x(i))^(deg-1);
        else
            val(i) = -deg*(2*e-x(i))^(deg-1);
        end
    end
end

function val = getf_der1(x, e, para, flag_type)
    if strcmp(flag_type, 'cos')
        val = getf_der1_cos(x, e, para);
    elseif strcmp(flag_type, 'power')
        val = getf_der1_power(x, e, para);
    end
end

function val = getf_der2_cos(x, a, para)
    N = length(x);
    val = zeros(1,N);
    for i =1:N
        if x(i)<=a && x(i)>=-a
            val(i) = -para^2* cos(para*x(i));
        elseif x(i)<-a
            val(i) = -para^2* cos(para*(-2*a-x(i)));
        else
            val(i) = -para^2* cos(para*(2*a-x(i)));
        end
    end
end

function val = getf_der2_power(x, e, deg)
    N = length(x);
    val = zeros(1,N);
    for i =1:N
        if x(i)<=e && x(i)>=-e
            val(i) = deg*(deg-1)*x(i)^(deg-2);
        elseif x(i)<-e
            val(i) = deg*(deg-1)*(-2*e-x(i))^(deg-2);
        else
            val(i) = deg*(deg-1)*(2*e-x(i))^(deg-2);
        end
    end
end

function val = getf_der2(x, e, para, flag_type)
    if strcmp(flag_type, 'cos')
        val = getf_der2_cos(x, e, para);
    elseif strcmp(flag_type, 'power')
        val = getf_der2_power(x, e, para);
    end
end


function val = getInt(coef_fh, e,right, M)
    coef_fh_exclude_first =  coef_fh(2:end);
    K = (1:M-1);
    val = (2*right/pi)*sum(coef_fh_exclude_first.*sin(K*pi*e/right)./K) + 2*e*coef_fh(1);
end

function val = mysimpson(x,y)
    y_middle = y(2:(end-1));
    y_middle_even = y_middle(2:2:end);
    y_middle_odd = y_middle(1:2:end);
    h = x(2)-x(1);
    val = (y(1)+y(end)+ 4*sum(y_middle_odd)+2*sum(y_middle_even))*h/3; %simpson 
end



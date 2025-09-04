%generate info for Table 2 in section 6 of the following papers
% On Trigonometric Interpolation and Its Applications I
% On Trigonometric Interpolation and Its Applications II

function main_sect6_table2()
    %
    % model = 'power', 'cos_alpha', 'sin_alpha', 'ln','abs', 'exp', 'cos'
    %
    model = 'power';
    m_powers = [6,8,10];
    num_for_plot = 2^12;   
    degs = [1,2];
    fig_type = '.fig';
    plot_location = 'northwest'; 
    do_graph = false;
    counter = 1;
    paras = degs;
    b = pi;
    if strcmp(model, 'sin_alpha') || strcmp(model, 'cos_alpha')
        paras = 4;
    end

    output = zeros(length(m_powers)*length(degs),9);    
    for para = paras
        for m_power = m_powers
            term = 2^m_power;
            M = term;
            N = 2*M;
            x = (0:(N-1))*2*pi/N-pi;
            x_plot = (0:(2*num_for_plot-1))*2*pi/(2*num_for_plot)-pi;
            fun = functor_fun(model, para, b);
            y = fun(x);
            [coef, delta] = cos_approx_engine_value2coef_wo_loop(y);
            filename = strcat('output/cos_approx_', model, '_M_', num2str(M), '_para_', num2str(para), fig_type);
   
            
            coef_standard = get_standard_coef(fun, M, pi);
            if do_graph
                [bound_1, bound_2] = get_der_estimation(model, para, b, M);
                plot_number = min(16, M);
                coef_estimation_x = 1:M;
                coef_estimation_x = coef_estimation_x(max(2,M-plot_number+2):end);
                title = get_title(model, para, M);
                legend_y = 'coef';
                legend_w = '$bounary$';
                xlabel = '$n$';
                ylabel = '$a_n$';
                filename_estimation = strcat('output/cos_approx_a_estimation_', model, '_der1_M_', num2str(M), '_para_', num2str(para), fig_type);              
                plot_latex_2(filename_estimation, coef_estimation_x, coef(coef_estimation_x), bound_2(coef_estimation_x), title, legend_y, legend_w, xlabel, ylabel, plot_location);
            end
            
            
            % diff_coef_cos_standard = max(abs(coef-coef_standard(1:M)));
            title = get_title(model, para, M);
            abs_cos = cos_approx_engine_plot(fun, coef, pi, num_for_plot, filename, title, title, title, title, plot_location, do_graph);
            abs_standard = cos_approx_engine_plot(fun, coef_standard, pi, num_for_plot, filename, title, title, title, title, plot_location, do_graph);
            funtor_der1 = functor_der_order1(model, para, b);
            der1_true = funtor_der1(x_plot);
            der1_cos = functor_der_order1_cos(coef, x_plot, b);
            der1_standard = functor_der_order1_cos(coef_standard, x_plot, b);
            funtor_der2 = functor_der_order2(model, para, b);
            der2_true = funtor_der2(x_plot);
            der2_cos = functor_der_order2_cos(coef, x_plot, b);
            der2_standard = functor_der_order2_cos(coef_standard, x_plot, b);
            if do_graph
                legend_z = '$y''_{M}(x)$';
                legend_y = '$y''(x)$';
                legend_w = '$y''_{std}(x)$';
                xlabel = '$x$';
                ylabel = '$y$';
                filename_der1 = strcat('output/main_sect6_table2_', model, '_der1_M_', num2str(M), '_para_', num2str(para), fig_type);
                title_d1 = strcat('The first derivative of ', '\quad ' , title);
                plot_latex_2(filename_der1, x_plot, der1_true, der1_cos, title_d1, legend_y, legend_z, xlabel, ylabel, plot_location);
                filename_der2 = strcat('output/main_sect6_table2_', model, '_der2_M_', num2str(M), '_para_', num2str(para), fig_type);
                legend_z = '$y"_{M}(x)$';
                legend_y = '$y"(x)$';
                legend_w = '$y"_{std}(x)$';                
                title_d2 = strcat('The second derivative of ', '\quad' ,title);
                plot_latex_2(filename_der2, x_plot, der2_true, der2_cos, title_d2, legend_y, legend_z, xlabel, ylabel, plot_location)
            end
            diff_der_1_cos_true = max(abs(der1_true-der1_cos));
            diff_der_1_standard_true = max(abs(der1_true-der1_standard));
            
            
            diff_der_2_cos_true = max(abs(der2_true-der2_cos));
            diff_der_2_standard_true = max(abs(der2_true-der2_standard));            
            %{
                max_cos_der_1 = max(abs(der1_cos));
                max_true_der_1 = max(abs(der1_true));
                max_standard_der_1 = max(abs(der1_standard));
                max_cos_der_2 = max(abs(der2_cos));
                max_true_der_2 = max(abs(der2_true));
                max_standard_der_2 = max(abs(der2_standard));
            %}
            output(counter,:)=[para, M, delta, abs_cos, abs_standard,   diff_der_1_cos_true, diff_der_1_standard_true, diff_der_2_cos_true, diff_der_2_standard_true];
            counter = counter+1;
            if do_graph
                if strcmp(model, 'power')
                    coef_N = real(ifft(y));
                    filename_N = strcat('output/main_sect6_table2_power_N_', num2str(N), '_deg_', num2str(para), fig_type);
                    alt = ones(1,N);
                    alt(2:2:N)= -1;
                    coef_N = coef_N.*alt;
                    cos_approx_engine_plot(fun, coef_N, pi, N, filename_N,title,title,title,title, plot_location, do_graph);
                    filename_2N = strcat('output/cos_approx_power_2N_', num2str(N), '_deg_', num2str(para), fig_type);
                    cos_approx_engine_plot(fun, coef_N, pi, 2*N, filename_2N,title,title,title,title, plot_location, do_graph);
                end
            end
        end
    end
    cos_approx_test_cos_file = strcat('output/sect6_tab2_', model, '.xlsx');
    vnames = {'deg','M', 'delta', 'err_fun', 'err_fun_standard',  'err_der1', 'err_der1_standard', 'err_der2', 'err_der2_standard'};
    mylib_writearray(vnames, output, cos_approx_test_cos_file)
end



function val =functor_der_order1_cos(coef, x, b)
    M = length(coef);
    order = (0:(M-1));
    n = length(x);
    val = zeros(1, n);
    for i = 1:n
        val(i)=-dot(coef.*order*(pi/b), sin(order*pi*x(i)/b));
    end
end

function val =functor_der_order2_cos(coef, x, b)
    M = length(coef);
    order = 0:(M-1);
    n = length(x);
    val = zeros(1, n);
    for i = 1:n
        val(i)=-dot(coef.*order.^2*(pi/b)^2, cos(order*pi*x(i)/b));
    end
end


function [title, title_d1, title_d2] = get_title(model, para, M)
     if strcmp(model, 'power')
         if para==1
            %title = strcat('$y(x)=1-(x/\pi)^2',' ,\quad M=', num2str(M), '$');
            title = strcat('$y(x)=1-(x/\pi)^2$');
         else
             %title = strcat('$y(x)=(1-(x/\pi)^2)^', num2str(para), ',\quad M=', num2str(M), '$');
             title = strcat('$y(x)=(1-(x/\pi)^2)^', num2str(para), '$');
         end
    elseif strcmp(model, 'cos_alpha')
        title = strcat('$y(x)=cos(', num2str(para), 'x), \quad M=', num2str(M), '$'); %para=alpha
    elseif strcmp(model, 'sin_alpha')
        title = strcat('$y(x)=sin(', num2str(para), 'x), \quad M=', num2str(M), '$');
    elseif strcmp(model, 'ln')
        if para == 1
            title = strcat('$y(x)=log(x^2+1)/log(\pi^2+1), \quad M=', num2str(M), '$');
        else
            title = strcat('$y(x)=(log(x^2+1))^', num2str(para), '/(log(pi^2+1))^', num2str(para), '\quad M=', num2str(M), '$');
        end
    elseif strcmp(model, 'abs')
        if para == 1
            title = strcat('$y(x)=(\pi-abs(x)), \quad M=', num2str(M), '$');
        else
            title = strcat('$y(x)=(\pi-abs(x))^', num2str(para), '\quad M=', num2str(M), '$');
        end        
    elseif strcmp(model, 'exp')
        if para == 1
            title = strcat('$y(x)=exp((\pi-abs(x)))/exp(\pi), \quad M=', num2str(M), '$');
        elseif para == 2
            title = strcat('$y(x)=exp((\pi^2-x^2))/exp(\pi^2), \quad M=', num2str(M), '$');
        else
            fprintf('invalid para %i', para)
        end        
    elseif strcmp(model, 'cos')
        title = strcat('$y(x)=\sum_{0\le j<', num2str(para), '}a_j cos(jx)', '\quad M=', num2str(M), '$');
    else
        fprintf('invalid model %s', model)
    end  
end

function [bound_1, bound_2] = get_der_estimation(model, para, b, M)
    N = 2*M;
    ll = 1:(M-1);
    cot_array = cot(ll*pi/N);
    x_N = -b + 0:(N-1)*(2*b)/N;
    fun_2 =functor_der_order2(model, para, b);
    D_2 = max(abs(fun_2(x_N)));
    fun_1 =functor_der_order1(model, para, b);
    D_1 = max(abs(fun_1(x_N)));
    bound_1 = (1+ cot_array)*(b*D_1*2^0.5/N);
    bound_1 = [0,bound_1];
    bound_2 = (1+ cot_array).^2*b^2*D_2*2*2^0.5/N^2;  
    bound_2 = [0,bound_2];
end



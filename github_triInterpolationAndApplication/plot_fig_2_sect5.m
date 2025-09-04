%generate Fig 2 and Fig 3 in section 5

function  plot_fig_2_sect5()
    fig_type = '.fig';
    plot_location = 'northwest'; 
    plot_size = 2^12;
    fun = @(x) (x-2.5).^2;
    delta = 1;
    left = 2;
    right = 3; 
    %q = 10;
    %qs = [4,5,6,7,8,9,10];
    qs = [10];
    color_code = 'blue';
    output = zeros(length(qs),2);
    counter = 1;
    for q = qs
        coef_yh = fourier_normalizer_cos_expansion_coef(fun, delta, left, right, q);
        D = right-left;
        x_plot = left-delta + (2*delta+D)/plot_size*(0:plot_size);  
        val = fourier_normalizer_coef2value_hf(coef_yh, left, right, delta, x_plot);
        
        y = fun(x_plot);
        %plot(x_plot, y, x_plot, val)
        x_range = x_plot(x_plot<= right & x_plot>=left);
        y_range = fun(x_range);
        val_range = fourier_normalizer_coef2value_hf(coef_yh, left, right, delta, x_range);
        output(counter,:) = [q, max(abs(y_range-val_range))];
        counter = counter + 1;
        if q == 10
            x_plot_period = left-delta-(D+2*delta) + 2*(2*delta+D)/plot_size*(0:plot_size);
            val_period = fourier_normalizer_coef2value_hf(coef_yh, left, right, delta, x_plot_period);
            filename = ['output/fig_3_sect5', num2str(q), fig_type];
            title = strcat('$f(x)=(x-2.5)^2$ vs $\hat{f}_M(x)$ over $[s-\delta,e+\delta]$ $(s,e,\delta)=(2,3,1)$');
            legend_y = '$f$';
            legend_z = '$\hat{f}_M$';
            xlabel_this = 'x';
            ylabel_this = 'y';
            plot_latex_2(filename, x_plot, y, val, title, legend_y, legend_z, xlabel_this, ylabel_this, plot_location )  
            file_name_period = ['output/fig_2_sect5_', num2str(q), fig_type];
            title_period = strcat('extend $f(x)$ to a period $[2s-e-3\delta,e+\delta]$,$(s,e,\delta)=(2,3,1)$');
            legend_period  = '';
            ylabel_period = 'y';
            plot_latex(file_name_period, x_plot_period, val_period, title_period, legend_period, xlabel_this, ylabel_period, plot_location, color_code)
        end
    end
    fourier_normalizer_test_hf_file = 'output/fig_2_sect5.xlsx';
    vnames = {'q', 'max_difference'};
    mylib_writearray(vnames, output, fourier_normalizer_test_hf_file)
end


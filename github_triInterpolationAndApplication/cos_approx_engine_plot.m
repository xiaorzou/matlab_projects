%plot function with cos expansion coef,  defined in [-b,b]

function error_abs = cos_approx_engine_plot(f, coef, b, N, filename, title_left, title_right, title_middle, title_full, plot_location, do_graph)
    addpath('D:/matlab/fourier_normalizer/fourier_normalizer_engine')  
    addpath('D:/matlab/cos_approx/cos_approx_engine')  
    addpath('D:/matlab/mylib')  
    term = length(coef);
    x = ((0:(N-1))*2*b/N-b);
    y = zeros(1,N);
    z = zeros(1,N);
    for k = 1:N
        y(k) = dot(coef, cos(((0:(term-1)).*x(k)*pi/b)));
        z(k) = f(x(k));
    end
    x_left = x(1:N/8);
    x_right = x((N/2-N/8+1):N/2);
    x_middle = x((N/4-N/16+1):(N/4+N/16));
    z_left = z(1:N/8);
    z_right = z((N/2-N/8+1):N/2);
    z_middle = z((N/4-N/16+1):(N/4+N/16));
    y_left = y(1:N/8);
    y_right = y((N/2-N/8+1):N/2);
    y_middle = y((N/4-N/16+1):(N/4+N/16));
    error_abs = max(abs(y-z));
    %title_left = 'f and cos approxmation over left wing';
    %title_right = 'f and cos approxmation over left right';
    %title_middle = 'f and cos approxmation over left middle';
    %title_full = 'f and cos approxmation';
    if do_graph
        legend_y = '$\bar{f}_N(x)$';
        legend_z = '$f(x)$';
        xlabel = '$x$';
        ylabel = '$y$';
        filename_full = strrep(filename,'.fig','_full.fig');
        filename_left = strrep(filename,'.fig','_left.fig');
        filename_right = strrep(filename,'.fig','_right.fig');
        filename_middle = strrep(filename,'.fig','_middle.fig');
        plot_latex_2(filename_left, x_left, y_left, z_left, title_left, legend_y, legend_z, xlabel, ylabel, plot_location)
        plot_latex_2(filename_right, x_right, y_right, z_right, title_right, legend_y, legend_z, xlabel, ylabel,plot_location)
        plot_latex_2(filename_middle, x_middle, y_middle, z_middle, title_middle, legend_y, legend_z, xlabel, ylabel,plot_location)
        plot_latex_2(filename_full, x, y, z, title_full, legend_y, legend_z, xlabel, ylabel,plot_location)
    end
end




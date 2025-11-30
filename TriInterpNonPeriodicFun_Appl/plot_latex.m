function plot_latex(file_name, X, Y, title_this, legend_this, xlabel_this, ylabel_this, location, color_code)
    close all
    low_y = min(Y);
    up_y = max(Y);
    low = low_y;
    up = up_y;
    xlim([X(1) X(end)])
    ylim([low up])
    %color_code = 'b',
    plot(X,Y,'-', 'Color', color_code)
    title1 = title(title_this,'Interpreter','latex');
    leg1 = legend({legend_this},'Interpreter','latex', 'Location', location); % R2018a and earlier
    xlab1 = xlabel(xlabel_this,'Interpreter','latex');
    ylab1 = ylabel(ylabel_this,'Interpreter','latex');
    set(leg1,'Interpreter','latex');
    set(xlab1,'Interpreter','latex');
    set(ylab1,'Interpreter','latex');
    set(title1,'Interpreter','latex');
    saveas(gcf,file_name);
end

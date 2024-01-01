function plot_latex_2(file_name, X, Y, Z, title_this, legend_y, legend_z, xlabel_this, ylabel_this )
    close all
    low_y = min(Y);
    low_z = min(Z);
    up_y = max(Y);
    up_z = max(Z);
    low = min([low_y, low_z]);
    up = max([up_y, up_z]);
    xlim([X(1) X(end)])
    ylim([low up])
    plot(X,Y,'--', X,Z,':')
    title1 = title(title_this,'Interpreter','latex');
    leg1 = legend({legend_y, legend_z},'Interpreter','latex'); % R2018a and earlier
    xlab1 = xlabel(xlabel_this,'Interpreter','latex');
    ylab1 = ylabel(ylabel_this,'Interpreter','latex');
    set(leg1,'Interpreter','latex');
    set(xlab1,'Interpreter','latex');
    set(ylab1,'Interpreter','latex');
    set(title1,'Interpreter','latex');
    saveas(gcf,file_name)
end


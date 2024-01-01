function plot_latex_3(file_name, X, Y, Z,W, title_this, legend_y, legend_z,legend_w, xlabel_this, ylabel_this )
    close all
    low_y = min(Y);
    low_z = min(Z);
    up_y = max(Y);
    up_z = max(Z);
    low_w = min(W);
    up_w = max(W); 
    low = min([low_y, low_z, low_w]);
    up = max([up_y, up_z, up_w]);
    axis([X(1) X(end) low up])
    xlim([X(1) X(end)])
    ylim([low up])
    plot(X,Y,'r', X,Z,'--', X,W,':')
    title1 = title(title_this,'Interpreter','latex');
    leg1 = legend({legend_y, legend_z, legend_w},'Interpreter','latex'); % R2018a and earlier
    xlab1 = xlabel(xlabel_this,'Interpreter','latex');
    ylab1 = ylabel(ylabel_this,'Interpreter','latex');
    set(leg1,'Interpreter','latex');
    set(xlab1,'Interpreter','latex');
    set(ylab1,'Interpreter','latex');
    set(title1,'Interpreter','latex');
    saveas(gcf,file_name)
end


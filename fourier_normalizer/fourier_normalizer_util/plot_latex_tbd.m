function plot_latex(file_name, X, Y, title_this, legend_this, xlabel_this, ylabel_this )
    close all
    plot(X,Y)
    if ~ strcmp(title_this,'')
        title1 = title(title_this,'Interpreter','latex');
    end
    if ~ strcmp(legend_this,'')
        leg1 = legend({legend_this},'Interpreter','latex'); % R2018a and earlier
        set(leg1,'Interpreter','latex');
    end
    xlab1 = xlabel(xlabel_this,'Interpreter','latex');
    ylab1 = ylabel(ylabel_this,'Interpreter','latex');
    set(xlab1,'Interpreter','latex');
    set(ylab1,'Interpreter','latex');
    set(title1,'Interpreter','latex');
    saveas(gcf,file_name)
end


function plot_latex_3(file_name, X, Y, Z,W, title_this, legend_y, legend_z,legend_w, xlabel_this, ylabel_this, location )
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
    plot(X,Y, '-' ,'Color', 	'r')
    hold on
    plot( X,Z,'-','Color',	'b')
    hold on 
    plot( X,W,'-' ,'Color',	'g')
    
    title1 = title(title_this,'Interpreter','latex');
    leg1 = legend({legend_y, legend_z, legend_w},'Interpreter','latex', 'Location', location); % R2018a and earlier
    xlab1 = xlabel(xlabel_this,'Interpreter','latex');
    ylab1 = ylabel(ylabel_this,'Interpreter','latex');
    set(leg1,'Interpreter','latex', 'FontName', 'Georgia');
    set(xlab1,'Interpreter','latex', 'FontName', 'Times');
    set(ylab1,'Interpreter','latex', 'FontName', 'Times');
    set(title1,'Interpreter','latex', 'FontName', 'Georgia');
    saveas(gcf,file_name)
end


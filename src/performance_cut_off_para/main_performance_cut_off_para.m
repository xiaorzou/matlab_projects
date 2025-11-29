
% the function generate table 1 of following paper
% "Trigonometric Interpolation on Non-Periodic Functions and its Applications"
%

function main_performance_cut_off_para()
    check_bdy_deg = 5;
    plot_size = 2^12;
    a = 1; 
    e = a;
    s = -e;
    isplot = false;
    test_type = 'default'; %reported in paper
    fig_type = '.fig';
    plot_location = 'northwest'; 
    if strcmp(test_type,'default')
        bs = [a+0.5*a];
        ds = [0.1,0.5,1];
        qs = [7,8];
    else
        bs = [2];
        ds = [0.5];
        qs = [8];
    end    
    rows = length(bs)*length(ds);
    deravative_at_pi_deg_n = zeros(rows,check_bdy_deg+4);
    deravative_at_0_deg_n = zeros(rows,check_bdy_deg+4);
    counter = 1;
    enhancement_status = zeros(rows,6);
    for bid = 1:length(bs)
        b=bs(bid);
        right = b;
        left = -right;
        x = fourier_normalier_get_grid(plot_size, 0, b);        
        for q = qs
            p = q-1;
            M = 2^(p);
            coef_n_impact = zeros(M,1);
            Aones2=ones(1,M);
            Aones2(2:2:end)=-ones(1,M/2);
            values = {[],[],[]};
            d_impact = containers.Map(ds, values);
            for did = 1:length(ds)
                d = ds(did);
                cut_off_para = d;
                h = fourier_normalizer_cut_off(x, left, s,e,right,cut_off_para);
                [h_coef_e, MaxError] = fourier_normalizer_get_coef_cut_off(q,left,s,e,right,cut_off_para);
                coef_n_impact(:, 1) = h_coef_e';
                h_der = cos_approx_engine_coef2value_enhance(h_coef_e, x, b, 1);
                %h_impact(:, 1) = h';
                d_impact(d) = h';
                if isplot && b == 2 
                    h_der_file_name = ['output/test_on_r_impact_', num2str(b),  '_', strrep(num2str(d), '.','_'), '_', num2str(q), '_cut_off' , fig_type];
                    title_der = ['$h(x)$ with $(b,d,q)=(', num2str(b), ',' num2str(d), ',', num2str(q),')$'];
                    legend_h_der = ['$h^{(1)}(x)$'];
                    xlabel = '$x$';
                    ylabel_h_der = ['$h^{(1)}$'];
                    plot_latex(h_der_file_name, x, h_der, title_der, legend_h_der, xlabel, ylabel_h_der, plot_location)
                    title_h = ['$h(x)$ and $h^{(', num2str(n), ')}(x)$ with $(b,q,d)=(',num2str(b), ',' num2str(q), ',', num2str(d), ')$'];
                    legend_h = '$h(x)$';
                    legend_h_der = ['$h^{(1)}(x)$'];
                    ylabel_h = '';
                    h_file_name = ['output/h_bdq',  num2str(b),  '_', strrep(num2str(d), '.','_'), '_', num2str(q), '_dim2_cut_off' , fig_type];
                    plot_latex_2(h_file_name, x, h, h_der, title_h, legend_h, legend_h_der, xlabel, ylabel_h, plot_location )
                end
                pib = pi*(0:1:M-1)/b;
                pib = pib.*pib;
                pib_2 = pib;
                alternative  = -1;
                deravative_at_pi_deg_n(counter,1:4) = [a,b,d,q];
                deravative_at_0_deg_n(counter,1:4) = [a,b,d,q];
                enhancement_status(counter,:) = [s,e,b-a,d,q,MaxError];
                for i=1:check_bdy_deg
                    deravative_at_pi_deg_n(counter,4+i) = log10(abs(alternative*dot(h_coef_e,pib.*Aones2)));
                    deravative_at_0_deg_n(counter,4+i) = log10(abs(alternative*dot(h_coef_e,pib)));
                    pib = pib.*pib_2;
                    alternative = -alternative;
                end
                counter = counter + 1;    
            end
            if q == qs(1) && b==2
                file_name_coef_d_impact = ['output/fig1_b_',  num2str(b), '_q_', num2str(q), '_cut_off' , fig_type];
                title_d = ['the impact on cut-off parameter'];
                legend_y = ['$h(x;', num2str(ds(1)) ,')$'];
                legend_z = ['$h(x;', num2str(ds(2)) ,')$'];
                legend_w = ['$h(x;', num2str(ds(3)) ,')$'];
                xlabel_this = 'x';
                ylabel_this = '$h$';
                s_pos = 1;
                e_pos = plot_size/4;
                y = d_impact(ds(1));
                z = d_impact(ds(2));
                w = d_impact(ds(3));
                plot_latex_3(file_name_coef_d_impact, x(s_pos:e_pos), y(s_pos:e_pos), z(s_pos:e_pos),w(s_pos:e_pos), title_d, legend_y, legend_z,legend_w, xlabel_this, ylabel_this, plot_location )
            end
        end
    end
    if strcmp(test_type,'default')
        vnames = {'a', 'b', 'r', 'q'};
        for deg = 1:check_bdy_deg
            vnames = [vnames , {strcat('deg_e',num2str(2*deg))}];
        end
        vnames_status = {'s', 'e', 'delta', 'r','q','MaxError'};        
        %file_name_left = 'output/test_on_r_impact_left.xlsx';
        %file_name_middle = 'output/test_on_r_impact_middle.xlsx';
        file_name_status = 'output/tab1.xlsx';
        %mylib_writearray(vnames, deravative_at_pi_deg_n, file_name_left)
        %mylib_writearray(vnames, deravative_at_0_deg_n, file_name_middle)
        mylib_writearray(vnames_status, enhancement_status, file_name_status)
    end
end

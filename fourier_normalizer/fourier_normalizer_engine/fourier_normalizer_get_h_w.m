% the function to get h(x) in the doc.
function output_table = fourier_normalizer_get_h_w(a,b,c,d,n,p,q, opt_plot, fig_type)
    N = 2^(q);
    M = 2^(p);  
    x = fourier_normalier_get_grid(N, 0, b);
    w_n_der_value = fourier_normalizer_get_w_der(x,a,b,c,d,n); %w^(n): nth derivative of w
    w_n_der_value = w_n_der_value/max(w_n_der_value);
      
    if mod(n,2)==0
        flag = 'sin';
    else
        flag = 'cos';
    end 
    w_n_der_coef = fourier_normalizer_value2coef(w_n_der_value, flag); 
    w_n_der_coef = w_n_der_coef(1:M); %fourier coefficient of w^(n)
    w_coef = fourier_normalizer_w_coef_transfer_deg_n_to_0(w_n_der_coef, n, b);       
    h_coef = fourier_normalizer_w2h(w_coef, b);
    saling = sum(h_coef);
    h_coef = h_coef/saling;
    %h_value = M*real(ifft(h_coef.*Aones2)); 
    h_value = fourier_normalizer_coef2value(h_coef, 'cos');
    
    %w_value_appr = M*imag(ifft(w_coef.*Aones2));
    w_value_appr = fourier_normalizer_coef2value(w_coef, 'sin');
    w_value_appr = w_value_appr/saling;
    w_value_appr_max = max(abs(w_value_appr));
    
 
    x_prime = fourier_normalier_get_grid(M, 0, b);
    w_value_true = fourier_normalizer_get_w_der(x_prime,a,b,c,d,0); 
    w_value_true = w_value_appr_max*w_value_true/max(abs(w_value_true));
    w_value_true = transpose(w_value_true);
    x_prime = transpose(x_prime);
    h_value = transpose(h_value);
    w_value_appr = transpose(w_value_appr);
    h_coef = transpose(h_coef);
    if opt_plot
        title_h = ['$h(x),$  $(a,b,c,d,n,p,q)=(' num2str(a) ',' num2str(b) ',' num2str(c) ',' num2str(d) ',' num2str(n) ',' num2str(p), ',' num2str(q) ')$'];
        title_w = ['$w(x),$  $(a,b,c,d,n,p,q)=(' num2str(a) ',' num2str(b) ',' num2str(c) ',' num2str(d) ',' num2str(n) ',' num2str(p), ',' num2str(q) ')$'];
        legend_h = '';%'$h(x)$';
        legend_appr = '$w_{appr}(x)$';
        legend_true = '$w_{true}(x)$';
        xlabel = '$x$';
        ylabel_h = '$h$';
        ylabel_w = '$w$';
        h_file_name = ['output/h_a_b_c_d_n_p_q_' num2str(a) '_' num2str(b) '_' num2str(c) '_' num2str(d) '_' num2str(n) '_' num2str(p) '_' num2str(q)  fig_type];
        w_file_name = ['output/w_a_b_c_d_n_p_q_' num2str(a) '_' num2str(b) '_' num2str(c) '_' num2str(d) '_' num2str(n) '_' num2str(p) '_' num2str(q)  fig_type];
        plot_latex(h_file_name,x_prime, h_value, title_h, legend_h, xlabel, ylabel_h)
        plot_latex_2(w_file_name, x_prime, w_value_appr, w_value_true, title_w, legend_appr, legend_true, xlabel, ylabel_w)
        %plot_latex(w_file_name,x_prime, w_value_appr, title_w, legend_w, xlabel, ylabel_w)
    end
    
    output_table = {table(x_prime, h_value, h_coef, w_value_appr, w_value_true),saling};
end






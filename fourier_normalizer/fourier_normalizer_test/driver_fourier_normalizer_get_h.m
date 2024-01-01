% the function to get h(x) in the doc.
function val = driver_fourier_normalizer_get_h()
    %{
        q = 10;
        p = 9;
        a = 1;
        c = 0.1;
        d = 0.1;
        a= 2;
        b = 3;
        check_bdy_deg = 5;
        flag_simpson = false;
        fig_type = '.fig'
    %}
    %cs = [0.01,0.05,0.1,0.5,1];
    cs = [0.1];
    %ds = [0.01,0.05,0.1,0.5,1];
    ds = [0.1];
    n = 2; % degree of derivaive of h to start with, can not be larger than 5
    fig_type = '.fig';
    flag_plot = true;
    p = 9;
    q = 10;
    a = 2;
    %bs = [3,5,7,10,15];
    bs = [3];
    N = 2^(q);
    M = 2^(p);
    Aones=ones(1,N);
    Aones(2:2:end)=-ones(1,N/2);
    Aones2=ones(1,M);
    Aones2(2:2:end)=-ones(1,M/2);
    
    for bid = 1:length(bs)
        b=bs(bid);
        X = fourier_normalier_get_fft_grid(N, 0, b);
        for cid = 1:length(cs)
            c=cs(cid);
            for did = 1:length(ds)
                d = ds(did);
                W_der_value = fourier_normalizer_get_w_der(X,a,b,c,d,n); %w^(n): nth derivative of w
                W_der_value = W_der_value/max(W_der_value);
                if mod(n,2)==0
                    W_der_coef = Aones.*(2*imag(ifft(W_der_value)));
                else
                    W_der_coef = Aones.*(2*real(ifft(W_der_value)));
                end
                W_der_coef = W_der_coef(1:M); %fourier coefficient of w^(n)
                W_coef_by_der = fourier_normalizer_get_w_from_derivative(W_der_coef, n, b);
                
                
                H_coef_by_n = fourier_normalizer_w2h(W_coef_by_der, b);
                saling = sum(H_coef_by_n);
                H_coef_by_n = H_coef_by_n/saling;
                %c_0 = -sum(F_H_by_n(2:end).*Aones2(2:end));
                %F_H_by_n =[c_0, F_H_by_n(2:end)];
                H_value_by_n = M*real(ifft(H_coef_by_n.*Aones2)); 

                w_appr_by_der = M*imag(ifft(W_coef_by_der.*Aones2));
                w_appr_by_der = w_appr_by_der/saling;
                

                W_0 = fourier_normalizer_get_w_der(X,a,b,c,d,0);
                W_0 = W_0/max(W_0);
                W_0_coef = Aones.*(2*imag(ifft(W_0)));
                W_0_coef = W_0_coef(1:M);
                H_0_coef = fourier_normalizer_w2h(W_0_coef, b);
                %c_0 = -sum(H_0_coef(2:end).*Aones2(2:end));
                %H_0_coef = [c_0, H_0_coef(2:end)];
                scaling_0 = sum(H_0_coef);
                % sum(A_h)=1 to meet the constrain h(0)=1!
                H_0_coef = H_0_coef/scaling_0;
                W_0_coef = W_0_coef/scaling_0;
                H_0_value = M*real(ifft(H_0_coef.*Aones2)); 
                w_appr_0 = M*imag(ifft(W_0_coef.*Aones2));

                x_prime = fourier_normalier_get_grid(M, 0, b);
                w_true = fourier_normalizer_get_w_der(x_prime,a,b,c,d,0);
                w_true = w_true/max(w_true)*max(w_appr_0);

                diff_w_n = max(abs(w_appr_by_der - w_true));
                diff_w_0 = max(abs(w_appr_0 - w_true));
                diff_H = max(abs(H_value_by_n - H_0_value));
                val = [diff_w_n, diff_w_0, diff_H];
                if flag_plot
                    %H = M*real(ifft(A_h.*Aones2));   
                    x_prime = fourier_normalier_get_grid(M, 0, b);
                    title_h = ['$h(x)$ $(a,b,c,d)=(' num2str(a) ',' num2str(b) ',' num2str(c) ',' num2str(d) ')$'];
                    title_w = ['$w(x)$ $(a, b,c,d)=(8\pi, 16\pi, 1)$'];
                    legend_h = '';%'$h(x)$';
                    legend_w = '';%'$w(x)$';
                    xlabel = '$x$';
                    ylabel_h = '$h$';
                    ylabel_w = '$w$';
                    h_file_name = ['output/h_a_b_c_d_q_p_' num2str(a) '_' num2str(b) '_' num2str(c) '_' num2str(d) '_' num2str(q) '_' num2str(p)  fig_type];
                    w_file_name = ['output/w_a_b_c_d_q_p_' num2str(a) '_' num2str(b) '_' num2str(c) '_' num2str(d) '_' num2str(q) '_' num2str(p)  fig_type];
                    plot_latex(h_file_name,x_prime, H_value_by_n, title_h, legend_h, xlabel, ylabel_h)
                    plot_latex(w_file_name,x_prime, w_true, title_w, legend_w, xlabel, ylabel_w)
                end
            end
        end
    end
end






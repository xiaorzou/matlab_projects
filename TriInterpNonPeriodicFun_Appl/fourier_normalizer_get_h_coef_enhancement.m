%{
 get coefficients of h's fourier expansion up to M-1 order, M=2^p
input:
    @p,b: the standard parameters in fourier_normalizer project, notice p=q-1
    @n: the degree of w's derivative,  n<=5 in current implementation, n
    has to be old to consider even function only. 
    @error_tol: error threshold in search coef in cos_approx_engine_value2coef, example = 10^(-12);
    @iter_num:  max iterate loops in searching coef in cos_approx_engine_value2coef, example 100;
output:
    @h_coef: the coefficients of fourier expanstion of h(x) from 0 up to M-1 
%}

%function [h_coef, delta, counter] = fourier_normalizer_get_h_coef_enhancement(a,b,c,d,n,q, error_tol, iter_num)
%further enhance using cos_approx_engine_value2coef_wo_loop to replace cos_approx_engine_value2coef
function [h_coef, delta] = fourier_normalizer_get_h_coef_enhancement(a,b,c,d,n,q)
    addpath('D:/matlab/cos_approx/cos_approx_engine') 
    N = 2^(q);
    x = fourier_normalier_get_grid(N, 0, b);
    w_n = fourier_normalizer_get_w_der(x,a,b,c,d,n); %w^(n): nth derivative of w
    w_n = w_n/max(w_n);
    if mod(n,2)==0
        disp('n has to be old')
        return
    end
    %opt = 'fix_imag';
    %[w_n_coef, delta, counter] = cos_approx_engine_value2coef(w_n, opt, error_tol, iter_num);
    [w_n_coef, delta] = cos_approx_engine_value2coef_wo_loop(w_n);
    w_coef = fourier_normalizer_w_coef_transfer_deg_n_to_0(w_n_coef, n, b);
    h_coef = fourier_normalizer_w2h(w_coef, b);
    saling = sum(h_coef);
    h_coef = h_coef/saling;  %make sure h(0)=1
end

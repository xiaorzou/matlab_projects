%{
 get coefficients of h's fourier expansion up to M-1 order, M=2^p
input:
    @p,q,b: the standard parameters in fourier_normalizer project
    @n: the degree of w's derivative,  n<=5 in current implementation
output:
    @h_coef: the coefficients of fourier expanstion of h(x) from 0 up to M-1 
%}

function h_coef = fourier_normalizer_get_h_coef(a,b,c,d,n,p,q)
    N = 2^(q);
    M = 2^(p);
    x = fourier_normalier_get_grid(N, 0, b);
    w_n = fourier_normalizer_get_w_der(x,a,b,c,d,n); %w^(n): nth derivative of w
    w_n = w_n/max(w_n);
    if mod(n,2)==0
        flag = 'sin';
    else
        flag = 'cos';
    end
    w_n_coef = fourier_normalizer_value2coef(w_n, flag);
    w_n_coef = w_n_coef(1:M);
    w_coef = fourier_normalizer_w_coef_transfer_deg_n_to_0(w_n_coef, n, b);
    h_coef = fourier_normalizer_w2h(w_coef, b);
    saling = sum(h_coef);
    h_coef = h_coef/saling;
end

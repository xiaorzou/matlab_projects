% get w(x)'s fourier coefficents from the coefficients of its n-th derivative
% w^(n)(x). notice w(x) is old function and cofficients =
% [0,B_1,...,B_{M-1}],  w_n_fourier_coefs=[0,..coef_{M-1}]
% input:  
% @w_n_fourier_coefs: the coef of fourier series of w^{(n)}(x)
% @n: degree of dervative of w(x)
% @b: the parameter b of fourier_normalizer parameter
% output:
% @val: w(x)'s fourier coefficents

function val = fourier_normalizer_w_coef_transfer_deg_n_to_0(w_n_fourier_coefs, n, b)
    k = floor(n/2);
    M = length(w_n_fourier_coefs);
    K = (b/pi)./(1:1:M-1);
    K2n = K.^(n);
    val = (-1)^(k)*K2n.*w_n_fourier_coefs(2:end);
    val = [0,val];
end


%{
 get h values at grid points x_k=-b + 2b*k/M, k=0,1,2,...,M-1
input:
    @a,b,c,d,n,p,q: the standard parameters in fourier_normalizer project
    @n: the degree of w's derivative,  n<=5 in current implementation
    @h_coef:  the coefficients of fourier expansion of h(x), term 0,...,M-1 
output:
    @h_value: h(x_k) for x_k=-b + 2b*k/M, k=0,1,2,...,M-1
%}

function h_value = fourier_normalizer_get_h_value(a,b,c,d,n,p,q,h_coef)
    M = 2^(p);
    Aones2=ones(1,M);
    Aones2(2:2:end)=-ones(1,M/2);
    if h_coef == false
        h_coef = fourier_normalizer_get_h_coef(a,b,c,d,n,p,q);
    end
    h_value = M*real(ifft(h_coef.*Aones2)); 
end

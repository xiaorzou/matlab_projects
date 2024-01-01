%{
 get f values at grid points x_k=-b + 2b*k/M, k=0,1,2,...,M-1
input:
    @a,b,c,d,n,p,q: the standard parameters in fourier_normalizer project
    @n: the degree of w's derivative,  n<=5 in current implementation
    @h_coef:  the coefficients of fourier expansion of h(x), term 0,...,M-1 
output:
    @h_value: h(x_k) for x_k=-b + 2b*k/M, k=0,1,2,...,M-1
%}

function value = fourier_normalizer_coef2value(coef, opt)
    M = length(coef);
    Aones2=ones(1,M);
    Aones2(2:2:end)=-ones(1,M/2);
    if strcmp(opt,'cos')
        value = M*real(ifft(coef.*Aones2));
    elseif strcmp(opt,'sin')
        value = M*imag(ifft(coef.*Aones2));
    else
        print('wrong opt in fourier_normalizer_coef2value')
        exit(0)
    end
     
end

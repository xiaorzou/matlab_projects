% f(x) = \sum_{1<=k<N}a_k*cos(kx), 
% g(x) = \sum_{0<=k<N}b_k*cos(kx), 
% f(x)*g(x) = \sum_{0<=k<2*N}c_k*cos(kx)
% N = length(a)
%input  
%   @a: coefficient array of f as above,  example: 1:10
%   @b: coefficient array of g as above,  example: 1:5
%output 
%   @c:  coefficents of f*g

function c = myfft_get_prod_fourier_coef(a, b)
    N = length(a);
    M = length(b);
    if N<M
        a = [a, zero(1,M-N)];
        N = M;
    elseif M<N
        b = [b, zero(1,N-M)];
    end
    c = zeros(1,2*N);
    c(1) = 0.5*a(1)*b(1) + 0.5*dot(a, b);

    for t = 2:N
        c(t)= 0.5*(dot(a(1:t), fliplr(b(1:t))) );
    end

    for t = N+1:2*N-1
        c(t)=  0.5*dot(a(t-N+1:N), fliplr(b(t-N+1:N)));
    end

    for t = 2: N
        c(t)= c(t) + 0.5*(dot(b(t:N),a(1:N-t+1)) + dot(a(t:N),b(1:N-t+1)));
    end
end
    





    
function test_myfft_get_prod_fourier_coef()
    a=10:-1:1;
    b=10:-1:1;
    N = length(a);
    M = length(b);
    if N<M
        a = [a, zero(1,M-N)];
        N = M;
    elseif M<N
        b = [b, zero(1,N-M)];
    end
    x = (1:1:10)/10;
    f = zeros(1,10);
    g = zeros(1,10);
    h = zeros(1,10);
    c = myfft_get_prod_fourier_coef(a, b);
    for i=1:10
        f(i)=dot(a, cos((0:1:N-1)*x(i)));
        g(i)=dot(b, cos((0:1:N-1)*x(i)));
        h(i)=dot(c, cos((0:1:(2*N-1))*x(i)))-f(i)*g(i);
    end
    max(abs(h))
end


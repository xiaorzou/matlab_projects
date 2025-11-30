% get the value of h^(deg) at x, where 
%   h(x)=coef(1)+coef(2)cos(pi*x/b)+... coef(N)cos((N-1)*x/b) 
% 

function y = cos_approx_engine_coef2value_enhance(coef, x, b, deg) 
    term = length(coef);
    N = length(x);
    y = zeros(1,N);
    factor = ((0:(term-1))*pi/b).^deg;
    for k = 1:N
        if mod(deg,2)==0
            l = deg/2;
            y(k) = (-1)^l*dot(coef.*factor, cos(((0:(term-1)).*x(k)*pi/b)));
        else
            l = (deg-1)/2;
            y(k) = (-1)^(l+1)* dot(coef.*factor, sin(((0:(term-1)).*x(k)*pi/b)));
        end
    end
end




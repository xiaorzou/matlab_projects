%plot function with cos expansion coef,  defined in [-b,b]

function y = cos_approx_engine_coef2value(coef, x, b) 
    term = length(coef);
    N = length(x);
    y = zeros(1,N);
    for k = 1:N
        y(k) = dot(coef, cos(((0:(term-1)).*x(k)*pi/b)));
    end
end




%plot function with cos expansion coef,  defined in [-b,b]

function y = cos_approx_engine_coef2value_general(coef, x, b, flag_trig) 
    term = length(coef);
    N = length(x);
    y = zeros(1,N);
    if strcmp(flag_trig,'cos')
        for k = 1:N
            y(k) = dot(coef, cos(((0:(term-1)).*x(k)*pi/b)));
        end
    elseif strcmp(flag_trig,'sin')
        for k = 1:N
            y(k) = dot(coef, sin(((0:(term-1)).*x(k)*pi/b)));
        end
    else
        disp(['not implemented with task ', task])
        return;
    end
end




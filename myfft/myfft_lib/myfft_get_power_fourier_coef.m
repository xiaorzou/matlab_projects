%return int^pi_{-pi} x^n g(x,m)dx,  g(x,m) = cos(mx) if opt = cos,
%g(x,m)=sin(mx) if opt = sin
function val = myfft_get_power_fourier_coef(n, m, opt)
    if n == 0 
        if strcmp(opt,'cos')
            if m == 0
                val = 2*pi;
            else
                val = 0;
            end
        else
            val = 0;
        end
    else
        if strcmp(opt,'cos')
            if m == 0
                val = (pi^(n+1)-(-pi)^(n+1))/(n+1);
            else
                val = -n*get_power_coef(n-1, m, 'sin')/m;
            end
        else
            val = (-1)^(m+1)*(pi^n-(-pi)^n)/m + (-1)^m*n*get_power_coef(n-1, m, 'cos')/m;
        end
    end
end

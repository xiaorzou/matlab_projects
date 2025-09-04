%{
input:
    @y: {y_1,...,y_N} with two cases and assciated flag
        y_k = \sum_{1\le j\le N} coef_j*cos((j-1)*pi*x_k/b), flag='cos'  
        y_k = \sum_{j\le j\le N} coef_j*sin((j-1)*pi*x_k/b), flag='sin'
    @flag: 'cos' or 'sin' as described above
output:
    coef: = {coef_1,...,coef_2} as described above

warning:
    in application,  one should use the first half of the coef to represent y=f(x),
    explained in the model doc, we use M (rathen N) terms in approximation:
    f(x) = \sum_{1\le j\le M} coef_j*cos((j-1)*pi*x/b) or
    f(x) = \sum_{1\le j\le M} coef_j*sin((j-1)*pi*x/b)
%}

function coef = fourier_normalizer_value2coef(y, flag)
    N = length(y);
    Aones=ones(1,N);
    Aones(2:2:end)=-ones(1,N/2);
    if strcmp(flag,'sin')
        coef = Aones.*(2*imag(ifft(y)));
    elseif strcmp(flag,'cos')
        coef = Aones.*(2*real(ifft(y)));
        %coef(1) = coef(1)/2; %newly added on 2024/1/9, not correct  remove
        %on June 4, 2025
        y_M = y(1:2:N); %add on june 4, 2025
        coef(1) = mean(y_M); %add on june 4, 2025
    else
        print('wrong flag in fourier_normalizer_get_coef_by_value');
        exit(0)
    end
end


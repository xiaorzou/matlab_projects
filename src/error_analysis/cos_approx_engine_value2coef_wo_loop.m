% y:  an array [y_1,...,y_N]  the function values f(x_k) at points 
%       x_k=-b+(k-1)*2b/N, k=1,2,...N, function f(x) is supposed to be periodic even function with period 2b. 
%       the function return the coef of cos expansion
%       f(x) ~ coef(1)+coef(2)*cos(1*pi*x/b)+coef(3)*cos(2*pi*x/b)+...+coef(M)*cos((M-1)*pi*x/b)

% ouput:
%   coef:  the targed coef as described above
%   delta: the maxium error between original function value and
%   approximated values over the points x_1, x_3, x_5, ...,  notice that
%   there are M=N/2 terms in the cos expansion. 

function [coef, delta] = cos_approx_engine_value2coef_wo_loop(y)
    addpath('D:/matlab/fourier_normalizer/fourier_normalizer_engine')  
    N = length(y);
    M = N/2;
    alt = ones(1,M);
    alt(2:2:M)= -1;
    coef = fourier_normalizer_value2coef(y, 'cos');
    coef = coef(1:M);
    y_M = y(1:2:N);
    coef(1) = mean(y_M);
    complex_fft = M*ifft(alt.*coef);
    y_fft_real = real(complex_fft);
    delta = max(abs(y_fft_real-y_M));
end
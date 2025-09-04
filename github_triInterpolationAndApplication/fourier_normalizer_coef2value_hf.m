
%{
input: coef_yh is the cos-expansion coefficient for yh(x), defined over
[-b,b] where b =  right-left + 2*delta;  where yh(x) is the even extension
from [0,b], where hy(x) = h(x)*y(x-b+right+delta), h(x) is cut-off defined
over [0,b] and =1 over [delta, b-delta]

output:  hf over x_plot, x_plot should be over the range [left-delta, right+delta]

%}
function  val = fourier_normalizer_coef2value_hf(coef_yh, left, right, delta, x_plot)
    %addpath('D:/matlab/fourier_normalizer/fourier_normalizer_engine')  
    %addpath('D:/matlab/cos_approx/cos_approx_engine')  
    D = right-left;
    %x_plot = left-delta + (2*delta+D)/plot_size*(0:plot_size);
    b =  D + 2*delta;
    x_shift = x_plot  - (right+delta) + b;
    val  = cos_approx_engine_coef2value(coef_yh, x_shift, b);
end

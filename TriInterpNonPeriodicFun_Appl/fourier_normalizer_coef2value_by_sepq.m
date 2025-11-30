
%{
input: coef_yh is the cos-expansion coefficient for yh(x), get by function 
fourier_normalizer_get_coef_by_fun(fun,s,e,p,q)
[s,e]:  the interval fun is supposed to be considered
p:   2^p: number of intervals that [s,e] to be divided into. 
q:   2^q: number of grid points before even extension.  as such,  total grid points used in cos algoirhtm is 2*2^q
notice that [s,e] is extended to [s-delta, e+delta] with delta =
[(2^q-2^p)/2] * [(e-s)/2^p]


defined over
[-b,b] where b =  r-l + 2*delta;  where 
    lambda = (e-s)/2^p;
    m = (2^q-2^p)/2;
    delta = lambda*m;
yh(x) is the even extension
from [0,b], where hy(x) = h(x)*y(x-b+right+delta), h(x) is cut-off defined
over [0,b] and =1 over [delta, b-delta]

output:  hf over x_plot, x_plot should be over the range [s-delta, e+delta]

%}
function  val = fourier_normalizer_coef2value_by_sepq(coef_yh, s,e,p,q,x_plot) 
    M = 2^q;
    n = 2^p;
    lambda = (e-s)/n;
    m = (M-n)/2;
    delta = lambda*m;
    l = s - delta;
    r = e + delta;
    %D = right-left;
    D = e-s;
    %x_plot = left-delta + (2*delta+D)/plot_size*(0:plot_size);
    b =  D + 2*delta;
    %x_shift = x_plot  - (right+delta) + b;
    x_shift = x_plot  - r + b;
    val  = cos_approx_engine_coef2value(coef_yh, x_shift, b);
end

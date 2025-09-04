%{
input: 
fun:  target fuction, which should be smooth in the range [left-detal, right+delta]
delta, left, right:   determine the smooth range of target function.
[left, right] determine the desired range that function value need to be
kept, [left-delta, left] and [right, right+delta] are buffter to define
cut-off function h. 

return: coefficients for y*h over the range (-b, b) where
b=(right-left)+2*delta, the periodic function is defined on [-b,b], based
on the even function, which is defined by f(x+left-delta) for x over
(0, b)
%}
function coef_yh = fourier_normalizer_cos_expansion_coef(fun, delta, left, right, q)
    addpath('D:/matlab/fourier_normalizer/fourier_normalizer_engine')  
    addpath('D:/matlab/cos_approx/cos_approx_engine')  
    global normalizer_d normalizer_q
    fourier_normalizer_const()
    N = 2^(q);
    x = left-delta + (right-left + 2*delta)/N*(0:(N-1));
    y = fun(x);
    h = fourier_normalizer_get_h_value_general(left,right,delta,normalizer_d,normalizer_q,x);
    yh = h.*y;
    yh = [0, fliplr(yh), yh(2:end)]; %even extention to the range with length 4*b
    [coef_yh, ~] = cos_approx_engine_value2coef_wo_loop(yh); 
end

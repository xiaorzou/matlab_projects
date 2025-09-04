%{
input:
q: N=2^q, number of grid points,  M=N/2 number of coefficients 
left/right:   h=0 outside [left, right]
s/e:  h=1 inside [s,e]
cuf_off_para: >0, the parameter of the cut-off generate function 
output: 
coef:   cos expansion coefficients with M=2^q/2 terms,  h(x) = \sum_{0\le j
<M} coef_j cos(j*pi*x/right)
error_this: max absolute error  
%}

function [coef, error_this] = fourier_normalizer_get_coef_cut_off(q,left,s,e,right,cut_off_para)
    addpath('D:/matlab/fourier_normalizer/fourier_normalizer_engine')  
    addpath('D:/matlab/cos_approx/cos_approx_engine')   
    N = 2^q;
    x_fit = left + (right-left)/N*(0:N-1);
    y_fit = fourier_normalizer_cut_off(x_fit, left, s,e,right, cut_off_para);
    [coef, ~] = cos_approx_engine_value2coef_wo_loop(y_fit);
    y_app = cos_approx_engine_coef2value(coef, x_fit, right);
    error_this = max(abs(y_fit-y_app));
end
    
    
    
    
    
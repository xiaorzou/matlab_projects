%{
input:   d, q:  parameter of h 
        left, right, delta:  define the shape of h,  h=1 over [left,
        right],  0 outside [left-delta, right+delta]

output:
return h cut-off function values of cut-off function h on the points x, that should be defined over the range [left-delta,
right+delta]. 
%}

function h_value = fourier_normalizer_get_h_value_general(left,right,delta,d,q, x)
    D = (right-left);
    b = delta + D/2;
    a = D/2;
    c = 1; %do not change it!
    n = 1; %do not change it!
    [h_coef, ~] = fourier_normalizer_get_h_coef_enhancement(a,b,c,d,n,q);
    x = x - (right+delta)+b;  % right+delta should be shifted to b!
    h_value = cos_approx_engine_coef2value(h_coef, x, b);
end

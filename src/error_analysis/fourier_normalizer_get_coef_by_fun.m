%{
return a_j, 0<=j<2^q
fun(x)*h(x) = sum_{j}coef_co_j cos(j*pi*(x-o)/b)
h(x): cut-off defined by s,e,p,q
x_0 = s-delta
n = 2^p
M = 2^q
2m+n = M
lamda = (e-s)/n
delta = m*lambda
o = s-delta
b = e+delta-o=e-s+2*delta
%}
function coef_co = fourier_normalizer_get_coef_by_fun(fun,s,e,p,q)
    global cuf_off_para 
    fourier_normalizer_const()
    M = 2^q;
    n = 2^p;
    lambda = (e-s)/n;
    m = (M-n)/2;
    delta = lambda*m;
    l = s - delta;
    r = e + delta;
    x = l + (r-l)/M*(0:M-1);
    y = fun(x);
    co = fourier_normalizer_cut_off(x, l,s,e,r,cuf_off_para);
    coy = co.*y;
    coy = [0, fliplr(coy), coy(2:end)]; %even extention to the range with length 4*b
    [coef_co, ~] = cos_approx_engine_value2coef_wo_loop(coy);
end

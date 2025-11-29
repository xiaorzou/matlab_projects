% input
% [left, s, e, right]:  left, start, end, right points.  
% output: a function such that f=0 outside [left, right], 1 within [s,e];


function val = fourier_normalizer_cut_off_base(x,cut_off_para)
% g(x)=f(x)/(f(x)+f(1-x)),  where f(x) is cut_off_generate
val = fourier_normalizer_cut_off_generator(x,cut_off_para)./(fourier_normalizer_cut_off_generator(x,cut_off_para)+fourier_normalizer_cut_off_generator(1-x,cut_off_para));
end


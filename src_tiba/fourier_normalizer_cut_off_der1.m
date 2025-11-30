% input
% [left, s, e, right]:  left, start, end, right points.  
% output: a function such that f=0 outside [left, right], 1 within [s,e];

function val = fourier_normalizer_cut_off_der1(x, left, s,e,right,cut_off_para)
% g((x-left)/(x-s))g(right-x)/(right-e))
    %val =  cut_off_base((x-left)/(s-left),cut_off_para).*cut_off_base((right-x)/(right-e),cut_off_para);
    val = 1/(s-left)*fourier_normalizer_cut_off_base_der1((x-left)/(s-left),cut_off_para).*fourier_normalizer_cut_off_base((right-x)/(right-e),cut_off_para) - 1/(right-e)*fourier_normalizer_cut_off_base((x-left)/(s-left),cut_off_para).*fourier_normalizer_cut_off_base_der1((right-x)/(right-e),cut_off_para);
end


function val = fourier_normalizer_cut_off_generator_der1(x,cut_off_para)
% e^{-cut_off_para/x} if x>0,  otherwise 0
    val = zeros(1,length(x));
    pos_p = find(x>0);
    %val(pos_p) = exp(-cut_off_para./x(pos_p));
    val(pos_p) = (2*cut_off_para./x(pos_p).^3).*exp(-cut_off_para./x(pos_p).^2);
end

function val = fourier_normalizer_cut_off_base_der1(x,cut_off_para)
% g(x)=f(x)/(f(x)+f(1-x)),  where f(x) is cut_off_generate
%val = cut_off_generator(x,cut_off_para)./(cut_off_generator(x,cut_off_para)+cut_off_generator(1-x,cut_off_para));
val = (fourier_normalizer_cut_off_generator_der1(x,cut_off_para).*(fourier_normalizer_cut_off_generator(x,cut_off_para)+fourier_normalizer_cut_off_generator(1-x,cut_off_para))...
    - fourier_normalizer_cut_off_generator(x,cut_off_para).*(fourier_normalizer_cut_off_generator_der1(x,cut_off_para) -  fourier_normalizer_cut_off_generator_der1(1-x,cut_off_para)))...
    ./(fourier_normalizer_cut_off_generator(x,cut_off_para)+fourier_normalizer_cut_off_generator(1-x,cut_off_para)).^2;
end


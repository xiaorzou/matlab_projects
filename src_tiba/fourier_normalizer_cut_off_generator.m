% input
% [left, s, e, right]:  left, start, end, right points.  
% output: a function such that f=0 outside [left, right], 1 within [s,e];


function val = fourier_normalizer_cut_off_generator(x,cut_off_para)
% e^{-cut_off_para/x} if x>0,  otherwise 0
    val = zeros(1,length(x));
    pos_p = find(x>0);
    %val(pos_p) = exp(-cut_off_para./x(pos_p));
    val(pos_p) = exp(-cut_off_para./x(pos_p).^2);
end

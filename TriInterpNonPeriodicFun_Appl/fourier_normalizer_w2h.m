% get fourier coefficients of h from those of w such that
% notice that H(b)=0, so we need coef_h_1-coef_h_2+...=0
% w(x) = coef_w(2)sin(pi*x/b)+coef_w(2)sin(2*pi*x/b)+... coef_w(M)sin((M-1)pi*x/b)
% output:
% h(x) =coef_h(1)+coef_h(2)cos(pi*x/b)+coef_h(3)cos(2*pi*x/b)+... coef_h(M)cos((M-1)pi*x/b)
% remark:
% h(-b) = coef_h(1)-coef_h(2)+coef_h(3)+... -coef_h(M) should be 0


function coef_h = fourier_normalizer_w2h(coef_w,b)
    M = length(coef_w);
    Aones2=ones(1,M);
    Aones2(2:2:end)=-ones(1,M/2);
    K = (b/pi)./(1:1:M-1);
    coef_h = -K.*coef_w(2:end);
    coef_h_0 = -sum(coef_h.*Aones2(2:end));
    coef_h =[coef_h_0, coef_h];
end

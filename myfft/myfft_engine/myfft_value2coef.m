%not complete yet,  2023-12-28

function [c_real, c_imag] = myfft_value2coef(f, opt)
%MYFFT_VALUE2COEF Summary of this function goes here
%   f_k=sum_{1<=j<=N}c_j*cos(2*pi*(j-1)*(k-1)/N) if opt='cos'
%   f_k=sum_{1<=j<=N}c_j*sin(2*pi*(j-1)*(k-1)/N) if opt='sin'
% 
%{
algorithm
assuming f=[f_1,f_2,...,f_N], construct complex array w with 
real = [f_1,f_2...f_{N/2},f_{N/2+1},f_{N/2},...,f_2]=f_1 + f(2:N/2)
+f(N/2+1)+f(N/2:2)
imag = [0,f_N,...,f_{N/2+2}, 0, -f_{N/2+2}, ..., -f_N]

%}
    N=length(f);
    myreal = [f(1:N/2+1), fliplr(f(2:N/2))];
    %myimage = [0,fliplr(f(N/2+2:N)),0, -f(N/2+2:N)];
    w = complex(myreal, myimage);
    c_complex = N*fft(w);
    c_real = real(c_complex);
    c_imag = imag(c_complex);
end


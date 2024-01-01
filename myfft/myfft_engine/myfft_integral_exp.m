% 
% \int^D(j)_{left)exp(h*x)q(x)dx, h=1,2,..,deg,  where D(j) is the grid
% points such that D(1)=left, D(N/2)=right. 
% input:
%   @density_coef: the coef of cos expansion of density function, 
%       output of my_fft_option_density_coef,  there are N terms
%   @deg:  the max degree to be applied in exp(h*x) as above
%   @left/right:  define the effective range of the density function q(x)
% output: double array of size (deg, N/2)


function FFTExp = myfft_integral_exp(density_coef, deg, left, right)
    N = length(density_coef);
    Aones=ones(1,N);
    Aones(2:2:end)=-ones(1,N/2);
    K=0:1:N-1;
    FFTExp=zeros(deg,N); 
    bma=right-left;
    lambda=bma/(N/2);
    for h=1:deg
        F2=density_coef./(1+ (pi/(h*bma))^2*K.^2);
        FFTExp(h,:)=-exp(left*h)*sum(F2)/h;
        expY=exp((2*left-right)*h)*exp(h*lambda*K);
        F2=F2.*Aones;
        FFTExp(h,:)=FFTExp(h,:)+N*expY.*real(ifft(F2))/h;
        F2=F2.*K;
        FFTExp(h,:)=FFTExp(h,:)+N*(pi/bma)*expY.*imag(ifft(F2))/h^2;
    end
    FFTExp=FFTExp(:, N/2+1:N);
end


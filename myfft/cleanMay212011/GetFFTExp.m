function void =GetFFTExp(F)
global FFTExp Left Right N Lambda
%global Sigma Rate Delta D %for testing purpose
%global F %keep this for testing purpose, we do not need this late

deg=1;

FFTExp=zeros(deg,N); % \int^D(j)_{Left)exp(h*x)q(x)dx, h=1,2,..,deg
bma=Right-Left;
bmapi=bma/pi;
%Y=0:1:(N-1);
%Y=2*bma*Y/N+2*Left-Right; %same as Y*Lambda+2*Left-Right
Aones=ones(1,N);
Aones(2:2:end)=-ones(1,N/2);
K=0:1:N-1;
%FFT for exponential function \int^D(j)_{Left)exp(h*x)q(x)dx, h=1,2,..,Deg
for h=1:deg
    F2=F./(1+ (pi/(h*bma))^2*K.^2);
    FFTExp(h,:)=-exp(Left*h)*sum(F2)/h;
    expY=exp((2*Left-Right)*h)*exp(h*Lambda*K);
    F2=F2.*Aones;
    FFTExp(h,:)=FFTExp(h,:)+N*expY.*real(ifft(F2))/h;
    F2=F2.*K;
    FFTExp(h,:)=FFTExp(h,:)+N*(pi/bma)*expY.*imag(ifft(F2))/h^2;
end
FFTExp=FFTExp(:, N/2+1:N);

function void =GetFFTG(F)
global FFTF FFTG Left Right N Lambda Deg M FFTAD
global Sigma Rate Delta D %for testing purpose
%global F %keep this for testing purpose, we do not need this late

FFTF=zeros(M,N); % Derivative of q(x) up to order M-1
FFTG=zeros(Deg,N); % \int^D(j)_{Left)exp(h*x)q(x)dx, h=1,2,..,Deg

bma=Right-Left;
bmapi=bma/pi;
Y=0:1:(N-1);
Y=2*bma*Y/N+2*Left-Right; %same as Y*Lambda+2*Left-Right
Aones=ones(1,N);
Aones(2:2:end)=-ones(1,N/2);
K=0:1:N-1;

%FFT for derivatives 
F1=Aones.*F;
FFTF(1,:)=N*real(ifft(F1));  % order 0 derivative q value

for i=1:((M-1)/2)
    F1=F1.*K*pi/(Right-Left);  % 
    FFTF(2*i,:)=(-1)^i*N*imag(ifft(F1)); % order 2*i-1 derivative
    F1=F1.*K*pi/(Right-Left);  %  
    FFTF(2*i+1,:)=(-1)^i*N*real(ifft(F1)); % order 2*i derivative
end
FFTF=FFTF(:,N/2+1:N);

%FFT for exponential function \int^D(j)_{Left)exp(h*x)q(x)dx, h=1,2,..,Deg
for h=1:Deg
    F2=F./(1+ (pi/(h*bma))^2*K.^2);
    FFTG(h,:)=-exp(Left*h)*sum(F2)/h;
    expY=exp((2*Left-Right)*h)*exp(h*Lambda*K);
    F2=F2.*Aones;
    FFTG(h,:)=FFTG(h,:)+N*expY.*real(ifft(F2))/h;
    F2=F2.*K;
    FFTG(h,:)=FFTG(h,:)+N*(pi/bma)*expY.*imag(ifft(F2))/h^2;
end
FFTG=FFTG(:, N/2+1:N);

function [FFTF, FFTAD] = FFT(F, N, M_A, M_D, Left, Right)
bma=Right-Left;
bmapi=bma/pi;
Aones=ones(1,N);
Aones(2:2:end)=-ones(1,N/2);
K=0:1:N-1;

%F=DensityCoeficient(model,inputVar, Left,Right, N);

%if 1==0
M = M_D;
FFTF=zeros(M,N); % Derivative of q(x) up to order M-1, M has to be odd!
F1=Aones.*F;
FFTF(1,:)=N*real(ifft(F1));  % order 0 derivative q value
%temp=FFTF(1,N/2+1:N);
for i=1:((M-1)/2)
    F1=F1.*K*pi/(Right-Left);  % 
    FFTF(2*i,:)=(-1)^i*N*imag(ifft(F1)); % order 2*i-1 derivative
    F1=F1.*K*pi/(Right-Left);  %  
    FFTF(2*i+1,:)=(-1)^i*N*real(ifft(F1)); % order 2*i derivative
end
FFTF=FFTF(:,N/2+1:N);

Y=0:1:(N-1);
Y=2*bma*Y/N+2*Left-Right; %same as Y*Lambda+2*Left-Right
MA = M_A;
FFTAD = zeros(MA*2,N);
Im = zeros (MA, N);
Re = zeros (MA, N);
F1_new=F(2:N)./K(2:N); %CDF
F1_new=[0,F1_new];
F1_new =bmapi*Aones.*F1_new;
for i=1:MA
    Im(i,:) = N*imag(ifft(F1_new));
    F1_new(2:N) = bmapi*F1_new(2:N)./K(2:N);
    Re(i,:) = N*real(ifft(F1_new));
    F1_new(2:N) = bmapi*F1_new(2:N)./K(2:N);
end
factor = 1;
for i=1:MA
    factor = factor *(2*i-1);
    FFTAD(2*i-1,:) = (-1)^(i-1)*Im(i,:) + (F(1)/factor)* (Y-Left).^(i*2-1);
    factor_2 = 1;
    for m=1:(i-1)
        factor_2 = factor_2*(2*m-1); 
        FFTAD(2*i-1,:) = FFTAD(2*i-1,:) + (-1)^(i-m-1)/factor_2*(Y-Left).^(m*2-1)*Re(i-m, N/2+1);
        factor_2 = factor_2 * 2*m ;
    end   
    factor = factor * 2*i;
    FFTAD(2*i,:) = (-1)^(i)*Re(i,:) + (F(1)/factor)* (Y-Left).^(i*2);
    factor_2 = 1;
    for m=1:i
        if m ~= 1
            factor_2 = factor_2*(2*m-2);
        end
        FFTAD(2*i,:) = FFTAD(2*i,:) + (-1)^(i-m)/factor_2*(Y-Left).^(m*2-2)*Re(i-m+1, N/2+1);
        factor_2 = factor_2 *(2*m-1);
    end       
end
FFTAD=FFTAD(:,N/2+1:N);
FFTAD=max(FFTAD,0);

%{
deg=1;
FFTExp=zeros(deg,N); % \int^D(j)_{Left)exp(h*x)q(x)dx, h=1,2,..,deg
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
%}
end


function Exp_fft = FFTExp(h,F, N, Left, Right)
%\int^D(j)_{Left)exp(h*x)q(x)dx
% F: density conefficient on values D(1), D(2), ....,D(N/2)
% D=(0:1:(N/2-1))*Lambda+Left; %N/2 elements
% D = L, L+Lambda, ..., R-Lambda = L+ (N/2-1)*Lambda
bma=Right-Left;
K=0:1:N-1;
Aones=ones(1,N);
Aones(2:2:end)=-ones(1,N/2);
%F=DensityCoeficient(model, Left,Right, N);
Lambda=(Right-Left)/(N/2);
Exp_fft=zeros(1,N); % 
F2=F./(1+ (pi/(h*bma))^2*K.^2);
Exp_fft(1,:)=-exp(Left*h)*sum(F2)/h;
expY=exp((2*Left-Right)*h)*exp(h*Lambda*K);
F2=F2.*Aones;
Exp_fft(1,:)=Exp_fft(h,:)+N*expY.*real(ifft(F2))/h;
F2=F2.*K;
Exp_fft(1,:)=Exp_fft(h,:)+N*(pi/bma)*expY.*imag(ifft(F2))/h^2;
Exp_fft=Exp_fft(:, N/2+1:N);
end



%fftf(order, j): q^(order-1)(D(j)), order =1, 2, 3, M
%fftf(1,j)=q(D(j))

function fftf=GetFFTF(model)
global Left Right N M
%global FFTF FFTG Left Right N Lambda Deg M FFTAD
%global Sigma Rate Delta D %for testing purpose

F=DensityCoeficient(model, Left,Right, N); 

fftf=zeros(M,N); % Derivative of q(x) up to order M-1

bma=Right-Left;
bmapi=bma/pi;
%Y=0:1:(N-1);
%Y=2*bma*Y/N+2*Left-Right; %same as Y*Lambda+2*Left-Right
Aones=ones(1,N);
Aones(2:2:end)=-ones(1,N/2);
K=0:1:N-1;

%FFT for derivatives 
F1=Aones.*F;
fftf(1,:)=N*real(ifft(F1));  % order 0 derivative q value

for i=1:((M-1)/2)
    F1=F1.*K*pi/(Right-Left);  % 
    fftf(2*i,:)=(-1)^i*N*imag(ifft(F1)); % order 2*i-1 derivative
    F1=F1.*K*pi/(Right-Left);  %  
    fftf(2*i+1,:)=(-1)^i*N*real(ifft(F1)); % order 2*i derivative
end
fftf=fftf(:,N/2+1:N);

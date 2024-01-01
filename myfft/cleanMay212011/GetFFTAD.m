%ffta(index, j)'=ffta(index-1,j), j=1,2,3,4,5.
%ffta(1,j)=CDF at $D(j)$

function fftad=GetFFTAD(modelType)
global Left Right N
%tic;
F=DensityCoeficient(modelType, Left,Right, N); 
bma=Right-Left;
bmapi=bma/pi;
Y=0:1:(N-1);
Y=2*bma*Y/N+2*Left-Right; %same as Y*Lambda+2*Left-Right
Aones=ones(1,N);
Aones(2:2:end)=-ones(1,N/2);
K=0:1:N-1;

fftad=zeros(5,N);

%first order anti, i.e. CDF
F1=F(2:N)./K(2:N);
F1=[0,F1];
F1=F1.*Aones;
fftad(1,:)=F(1)*(Y-Left)+N*bma*imag(ifft(F1))/pi;
%fftad=fftad(N/2+1:N);

%second order anti
F1(2:N)=F1(2:N)./K(2:N);
temp2=N*real(ifft(F1));
fftad(2,:)=0.5*F(1)*(Y-Left).^2-bmapi^2.*temp2+ temp2(N/2+1)*bmapi^2;

%third order
F1(2:N)=F1(2:N)./K(2:N);
fftad(3,:)=(1/6)*F(1)*(Y-Left).^3-bmapi^3*(N*imag(ifft(F1)))+(Y-Left)*bmapi^2*temp2(N/2+1);

%fourth order
F1(2:N)=F1(2:N)./K(2:N);
temp3=N*real(ifft(F1));
fftad(4,:)=(F(1)/24)*(Y-Left).^4+bmapi^4*temp3-bmapi^4*temp3(N/2+1)+0.5*(Y-Left).^2*bmapi^2*temp2(N/2+1);

%fifth order
F1(2:N)=F1(2:N)./K(2:N);
temp4=N*imag(ifft(F1));
fftad(5,:)=(F(1)/120)*(Y-Left).^5+bmapi^5*temp4-bmapi^4*temp3(N/2+1)*(Y-Left)+bmapi^2*temp2(N/2+1)*(Y-Left).^3/6;

fftad=fftad(:,N/2+1:N);

%sprintf('in GetAntiDerivative, abs(CDF-1) is appxoimated by %f',
%abs(fftad(N/2)-1))
%t1=toc;
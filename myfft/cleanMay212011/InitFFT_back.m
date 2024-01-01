%   PIntA(deg, j)=\int^{D(j)+DeltaA}_{D(j)}(x-D(j))^(deg-1)q(x)dx
%   deg=1,2,3,...,Deg+1, j=1,2,...,N/2
%   \int^{D(N/2+1)}_{D(N/2))(x-D(j))^(deg-1)q(x)dx is set to be 0 for all deg=1,2,Deg+1

function InitFFT_back(model) 
global N Deg D Left Right Lambda  FFTExp ModelScaling FFTAD %FFTF
global SA SB
global Rate Sigma Delta  NumA NumB B GP CP D3

bma=Right-Left;
bmapi=bma/pi;
Aones=ones(1,N);
Aones(2:2:end)=-ones(1,N/2);
K=0:1:N-1;

F=DensityCoeficient(model, Left,Right, N);

if 1==0
M=1
FFTF=zeros(M,N); % Derivative of q(x) up to order M-1, M has to be odd!
F1=Aones.*F;
FFTF(1,:)=N*real(ifft(F1));  % order 0 derivative q value
temp=FFTF(1,N/2+1:N);
for i=1:((M-1)/2)
    F1=F1.*K*pi/(Right-Left);  % 
    FFTF(2*i,:)=(-1)^i*N*imag(ifft(F1)); % order 2*i-1 derivative
    F1=F1.*K*pi/(Right-Left);  %  
    FFTF(2*i+1,:)=(-1)^i*N*real(ifft(F1)); % order 2*i derivative
end
FFTF=FFTF(:,N/2+1:N);
end







%init FFTAD
Y=0:1:(N-1);
Y=2*bma*Y/N+2*Left-Right; %same as Y*Lambda+2*Left-Right
FFTAD=zeros(Deg+1,N);

F1=F(2:N)./K(2:N); %CDF
F1=[0,F1];
F1=F1.*Aones;
FFTAD(1,:)=F(1)*(Y-Left)+N*bma*imag(ifft(F1))/pi;

F1(2:N)=F1(2:N)./K(2:N); %second order anti
temp2=N*real(ifft(F1));
FFTAD(2,:)=0.5*F(1)*(Y-Left).^2-bmapi^2.*temp2+ temp2(N/2+1)*bmapi^2;

F1(2:N)=F1(2:N)./K(2:N); %third order
FFTAD(3,:)=(1/6)*F(1)*(Y-Left).^3-bmapi^3*(N*imag(ifft(F1)))+(Y-Left)*bmapi^2*temp2(N/2+1);

F1(2:N)=F1(2:N)./K(2:N); %fourth order
temp3=N*real(ifft(F1));
FFTAD(4,:)=(F(1)/24)*(Y-Left).^4+bmapi^4*temp3-bmapi^4*temp3(N/2+1)+0.5*(Y-Left).^2*bmapi^2*temp2(N/2+1);

F1(2:N)=F1(2:N)./K(2:N); %fifth order
temp4=N*imag(ifft(F1));
FFTAD(5,:)=(F(1)/120)*(Y-Left).^5+bmapi^5*temp4-bmapi^4*temp3(N/2+1)*(Y-Left)+bmapi^2*temp2(N/2+1)*(Y-Left).^3/6;

if Deg>4
    error('not implemented yet, in InitFFT');
end


FFTAD=FFTAD(:,N/2+1:N);
FFTAD=max(FFTAD,0);







%D3(base,deg,pos)=\int^(D(B(base+1))-D(pos))_{D(B(base))-D(pos)) (w+D(pos)-D(B(base))^(Deg-deg+1)q_W(w)dw
%base: 1:NumA+NumB
%deg: 1:Deg
%pos: 1:(NumA+NumB)*Deg+1

%tic;
D3=zeros(NumA+NumB,Deg,(NumA+NumB)*Deg+1);
for degg=1:Deg
    deg=Deg-degg+1;
    for base=1:(NumA+NumB)
        poslow=B(base)-GP+CP;
        poslow=min(max(poslow,1),N/2);
        posup=B(base+1)-GP+CP;
        posup=min(max(posup,1),N/2);
        %D3(base,degg,:)
        temp=FFTAD(1,posup).*(D(posup)+D(GP)-D(B(base))).^deg-FFTAD(1,poslow).*(D(poslow)+D(GP)-D(B(base))).^deg ;
        factor=1;
        for k=1:deg
            factor=factor*(deg-k+1);
            temp=temp+(-1)^k*factor*( FFTAD(k+1,posup).*(D(posup)+D(GP)-D(B(base))).^(deg-k)-...
                FFTAD(k+1,poslow).*(D(poslow)+D(GP)-D(B(base))).^(deg-k));
        end
        D3(base,degg,:)=max(temp,0);
    end
end

%t=toc;
%sprintf('time is %f,', t)

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

if 1==0
[FFTF(1),FFTF(N/2),1-FFTAD(1,N/2), MyChar(model,-complex(0,1)*1)-FFTExp(N/2)]
%FFTF'
exit(0)
end




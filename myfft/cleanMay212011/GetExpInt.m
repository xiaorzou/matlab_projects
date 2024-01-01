%fftg(h,j)=\int^(D(j))_{Left} e^{i*x)q(x)dx,  h=1,2,3,...,Deg, j=1,2,... N/2
function fftg=GetExpInt(model)
global Left Right N M Deg

F=DensityCoeficient(model, Left,Right, N); 
fftg=zeros(Deg,N);
bma=Right-Left;
bmapi=bma/pi;
%Y=0:1:(N-1);
%Y=2*bma*Y/N+2*Left-Right; %same as Y*Lambda+2*Left-Right
Aones=ones(1,N);
Aones(2:2:end)=-ones(1,N/2);
K=0:1:N-1;
%FFT for exponential function \int^D(j)_{Left)exp(h*x)q(x)dx, h=1,2,..,Deg
for h=1:Deg
    F2=F./(1+ (pi/(h*bma))^2*K.^2);
    fftg(h,:)=-exp(Left*h)*sum(F2)/h;
    expY=exp((2*Left-Right)*h)*exp(h*Lambda*K);
    F2=F2.*Aones;
    fftg(h,:)=fftg(h,:)+N*expY.*real(ifft(F2))/h;
    F2=F2.*K;
    fftg(h,:)=fftg(h,:)+N*(pi/bma)*expY.*imag(ifft(F2))/h^2;
end
fftg=fftg(:, N/2+1:N);
disp('in GetExpInt')
[abs(fftg(1,N/2)-MyChar(modelType,-complex(0,1)*1)),
 abs(fftg(2,N/2)-MyChar(modelType,-complex(0,1)*2)), 
 abs(fftg(3,N/2)-MyChar(modelType,-complex(0,1)*3)), 
 abs(fftg(4,N/2)-MyChar(modelType,-complex(0,1)*4))]
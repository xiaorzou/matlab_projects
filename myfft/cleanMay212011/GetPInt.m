%   PInt(deg, j)=\int^{D(j+1)}_{D(j)}(x-D(j))^(deg-1)q(x)dx
%   deg=1,2,3,...,Deg+1, j=1,2,...,N/2
%   \int^{D(N/2+1)}_{D(N/2))(x-D(j))^(deg-1)q(x)dx is set to be 0 for all deg=1,2,Deg+1

function pint=GetPInt(model) 
global N Deg D Left Right Lambda FFTF 

FFTF=GetDerivative(model);
pint=zeros(Deg+1,N/2-1);
%up=D(2:N/2);
%low=D(1:N/2-1);

ffta=GetAntiDerivative(model);
for deg=1:Deg+1
    pint(deg,:)=ffta(1,2:N/2)*Lambda^(deg-1)-ffta(1,1:N/2-1)*0^(deg-1);
    factor=1;
    for k=1:deg-1
        factor=factor*(deg-k) ;
        pint(deg,:)=pint(deg,:)+(-1)^k*factor*(ffta(k+1,2:N/2)*Lambda^(deg-k-1)-...
            ffta(k+1,1:(N/2-1))*0^(deg-k-1) );
    end
end
pint=[pint, zeros(Deg+1,1)];




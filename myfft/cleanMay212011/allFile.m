
%%%%%%%%%%%   from final fold %%%%%%%%%%%%%%%%%%%%%

%x do not have to be in grid
%S=exp(x): the price of stock
%val=continuation value
%di: (i-1)th derivative
%vector =[d1, d2,....d(deg)]

%if option ==1, find contination, first derivative and second derivative
%if option ==0, just continuation

function vector=Continuation(pos, option)
global KStar Oalpha Obeta Strike Boundary Discount NumFun Coef NumPeriods EpsilonSearching Lambda
global CP D Deg Euro
global Index Points KP Left Right  SStar %donot need this one late.

if x<Left | x>Right
    error('wrong range for x in contiaiton');
end
x=D(pos);
val=0;
s=exp(x);

low=min(max(KStar-x, Left), Right);
%up1=min(max(log(Strike)-x, Left), Right); %log(Strike)=0
%up2=min(max(D(Boundary(2))-x, Left), Right);

cdflow=GetDerivative(low,0); %CDF(low)
qlow=GetDerivative(low,1);   %q(low)
qdlow=GetDerivative(low,2);  %q'(low)
explow=GetExpF(low);

part1=Oalpha*Strike*cdflow-Obeta*s*explow;

d1=0;
d2=0;
part1=0;
part2=0;
delta=D(Boundary(2))-D(Boundary(1));
if option==1
    d1=(polyval(Coef(1,:),KStar)-Oalpha*Strike+Obeta*SStar)*qlow - Obeta*s*explow;
    coed=Coef(1,:).*(Deg:-1:0);
    coed=coed(1:Deg)/delta;
    d2=( Obeta*s + polyval(coed, KStar))*qlow-...    
    (polyval(Coef(1,:),KStar)-Oalpha*Strike+Obeta*SStar)*qdlow-...
    Obeta*s*explow;
end


firstPartCon=GetContP1Int(x); 
part2=sum(Coef(1,:).firstPartCon);
if option==1
    d1=d1+sum(coed.*firstPartCon(1:Deg));
    coedd=coed.*(Deg-1:-1:0)/delta;
    d2=d2+sum(coedd(1:Deg-1).*firstPartCon(1:Deg-1));
end

for j=2:NumFun
    conPart2=GetContP2Int(pos,j);
    part2=part2+sum(Coef(j,:).conPart2);
    if option==1
        delta=D(Boundary(j+1))-D(Boundary(j));
        coed=Coef(j,:).*(Deg:-1:0);
        coed=coed(j, 1:Deg)/delta;
        d1=d1+sum(coed.*conPart2(1:Deg));
        coedd=coed.*(Deg-1:-1:0)/delta;
        d2=d2+sum(coedd(1:Deg-1).*conPart2(1:Deg-1));
    end
end
val=Discount*(part1+part2);
vector=[max(val,0), d1, d2];
if val <0
    val
end

%[\int^{d2-x}_{kstar-x} [(u+x-d(1))/(d(2)-d(1))]^deg q(u)du, 
%\int^{d2-x}_{kstar-x} [(u+x-d(1))/(d(2)-d(1))]^(deg-1) q(u)du,...,
%\int^{d2-x}_{kstar-x} [(u+x-d(1))/(d(2)-d(1))]^0 q(x),]
function firstPartCon=GetContP1Int(x)
global Boundary D Fiba Deg

delta=D(Boundary(2))-D(Boundary(1));


function conPart2=GetContP2Int(pos,j);
global Boundary D Fiba Deg FFTAD

position=Boundary(j+1)-pos+N/4+1;
F=FFTAD(pos,1:Deg+1);
vect=-1;
for p=1:Deg+1
    position(Deg+2-p)=F(1);
    coe=1;
    for i=1:p
        coe=coe*(p-i+1);
        position(Deg+2-p)=position(Deg+2-p)+(-1)^i*F(1+i)*coe; 
    end
end
----------------------------

function coefdensity=DensityCoeficient(modelType, leftend,rightend, N)
x=0:1:(N-1);
coefdensity=MyChar(modelType, x*pi/(rightend-leftend)).*exp(-x*leftend*pi/(rightend-leftend)*complex(0,1));
coefdensity=2*real (coefdensity)/(rightend-leftend);
coefdensity(1)=coefdensity(1)/2;

---------------------------------

%q_{(order)}(at): if order=0, the value of q(at), otherwise the value of
%antiderivative at 'at', expend around D(lpos)
%order =0, 1, 2, ..., 5 
%q(at)=GetAntiDerivative(at, lpos, 0)
%\int^at_{Left} q(x)dx=GetAntiDerivative(at, lpos, 1)

function val=GetAntiDerivative(x, order)
global D FFTF FFTAD Lambda

if x<Left | x>Right
    error('Wrong range for x %d', x)
end
q=floor((x-Left)/Lambda);
r=x-Left - Lambda *q;
if q==N/2
    q=q-1;
end

lpos=q+1;
if r> (Lambda/2)
    lpos=lpos+1;
end

delta=x-D(lpos);

if abs(delta)>Lambda
    error('in GetAntiDerivative')
end
if order==0
    if abs(delta)<FloatingError
        val=FFTF(1, lpos);
    else
        val=FFTF(1, lpos)+delta*FFTF(2,lpos)+0.5*delta^2*FFTF(3,lpos)+delta^3*FFTF(4,lpos)/6+...
            delta^4*FFTF(5, lpos)/24 +delta^5*FFTF(6, lpos)/120 +delta^6*FFTF(7, lpos)/720;
    end
elseif order ==1 % \int^at_{-\infty} q(x)dx
    if abs(delta)<FloatingError
        val=FFTAD(1,lpos);
    else
        val=FFTAD(1,lpos)+delta*FFTF(1,lpos)+delta^2*FFTF(2,lpos)/2+delta^3*FFTF(3,lpos)/6+delta^4*FFTF(4, lpos)/24 +...
        delta^5*FFTF(5, lpos)/120+delta^6*FFTF(6, lpos)/720+delta^7*FFTF(7, lpos)/5040;
    end
elseif order==2
    if abs(delta)<FloatingError
        val=FFTAD(2,lpos);
    else
        val=FFTAD(2,lpos)+delta*FFTAF(1,lpos)+delta^2*FFTF(1,lpos)/2+delta^3*FFTF(2,lpos)/6 + ...
        delta^4*FFTF(3, lpos)/24 +delta^5*FFTF(4, lpos)/120 +...
        delta^6*FFTF(5, lpos)/720+delta^7*FFTF(6, lpos)/5040+delta^8*FFTF(7, lpos)/40320;
    end
elseif order==3
    if abs(delta)<FloatingError
        val=FFTAD(3,lpos);
    else
        val=FFTAD(3,lpos)+delta*FFTAD(2,lpos)+delta^2*FFTAD(1,lpos)/2+delta^3*FFTF(1,lpos)/6 + ...
        delta^4*FFTF(2, lpos)/24 +delta^5*FFTF(3, lpos)/120 +...
        delta^6*FFTF(4, lpos)/720+delta^7*FFTF(5, lpos)/5040+delta^8*FFTF(6, lpos)/40320+...
        delta^9*FFTF(7, lpos)/362880;
    end
elseif order==4
    if abs(delta)<FloatingError
        val=FFTAD(4,lpos);
    else
        val=FFTAD(4,lpos)+delta*FFTAD(3,lpos)+delta^2*FFTAD(2,lpos)/2+delta^3*FFTAD(1,lpos)/6 + ...
        delta^4*FFTF(1, lpos)/24 +delta^5*FFTF(2, lpos)/120 +...
        delta^6*FFTF(3, lpos)/720+delta^7*FFTF(4, lpos)/5040+delta^8*FFTF(5, lpos)/40320+...
        delta^9*FFTF(6, lpos)/362880+delta^10*FFTF(7, lpos)/3628800;
    end
elseif order==5
    if abs(delta)<FloatingError
        val=FFTAD(5,lpos);
    else
        val=FFTAD(5,lpos)+delta*FFTAD(4,lpos)+delta^2*FFTAD(3,lpos)/2+delta^3*FFTAD(2,lpos)/6 + ...
        delta^4*FFTAD(1, lpos)/24 +delta^5*FFTF(1, lpos)/120 +...
        delta^6*FFTF(2, lpos)/720+delta^7*FFTF(3, lpos)/5040+delta^8*FFTF(4, lpos)/40320+...
        delta^9*FFTF(5, lpos)/362880+ delta^10*FFTF(6, lpos)/3628800 +delta^11*FFTF(7, lpos)/39916800;
    end
else
    error('not implemented yet, anti derivative order up to 5 in GetAntiDerivative');
end
-----------------------------------------

%q_{(order)}(at): if order=0, the value of q(at), otherwise the value of
%antiderivative at 'at', expend around D(lpos)
%order =0, 1, 2, ..., 5 
%q(at)=GetAntiDerivative(at, 0)
%\int^at_{Left} q(x)dx=GetAntiDerivative(at,1)

function val=GetCDF(x)
global D FFTF CDF Lambda

if x<Left | x>Right
    error('Wrong range for x %d', x)
end
q=floor((x-Left)/Lambda);
r=x-Left - Lambda *q;
if q==N/2
    q=q-1;
end

lpos=q+1;
if r> (Lambda/2)
    lpos=lpos+1;
end

delta=x-D(lpos);

if abs(delta)>Lambda
    error('in GetAntiDerivative')
end

if abs(delta)<FloatingError
    val=CDF(1,lpos);
else
    val=CDF(1,lpos)+delta*FFTF(1,lpos)+delta^2*FFTF(2,lpos)/2+delta^3*FFTF(3,lpos)/6+delta^4*FFTF(4, lpos)/24 +...
        delta^5*FFTF(5, lpos)/120+delta^6*FFTF(6, lpos)/720+delta^7*FFTF(7, lpos)/5040;
end

------------------------------------
%  order=0, 1, 2, 3, ... M
%  GetDerivative(x, 0)=CDF(x)
%  GetDerivative(x, 1)=q(x)
%  GetDerivative(x, order)=GetDerivative(x, order-1)'
function val=GetDerivative(x, order)
global D FFTF Lambda M Left Right N FloatingError CDF

if order>M
    error('in GetDensity, the order %d is bigger than M:%d ', order, M);
end
if (x<Left | x>Right)
    if order==1
        val=0; %density
    elseif order==0 %CDF
        if x>Right 
            val=1;
        elseif x<Left
            val=0;
        end
    else %for other order, we are not sure how to set the value.
        error('wrong range for x in GetDerivative');
    end
    return;
end
q=floor((x-Left)/Lambda);
r=x-Left - Lambda *q;
if q>=N/2
    q=N/2-1;
    r=0;
end
if (r/Lambda)<FloatingError
%    disp('enter Gend1');
    if order==0
        val=CDF(q+1);
    else
        val=FFTF(order,q+1);
    end
    return;
elseif (r/Lambda)>(1-FloatingError)
    if order==0
        val=CDF(q+2);
    else
        val=FFTF(order, q+2);
    end
    return;
else
    pos=0;
    if (r/Lambda)<0.5
        pos=q+1;
    else
        pos=q+2;
    end
    delta=x-D(pos);
    c=1;
    if order==0 %for cdf
        val=CDF(pos);
        %+delta*FFTF(1,lpos)+delta^2*FFTF(2,lpos)/2+delta^3*FFTF(3,lpos)/6+delta^4*FFTF(4, lpos)/24 +...
        %delta^5*FFTF(5, lpos)/120+delta^6*FFTF(6, lpos)/720+delta^7*FFTF(7, lpos)/5040;
    else
        val=FFTF(order,pos);
    end
    for i=1:M-order
    %for i=1:6
        c=c*delta/i;
        val=val+c*FFTF(i+order,pos);
    end
end
---------------------------------
%int^x_Left e^(u) q(u)du   for x\in [Left Right]
function val = GetExpF(x)
global N Left Right Lambda  FFTF FFTExp D M  FloatingError
if x<Left | x>Right
    error('Wrong range for x %d', x)
end
q=floor((x-Left)/Lambda);
r=x-Left - Lambda *q;
if q==N/2
    q=q-1;
end
if abs(r/Lambda)<FloatingError
    val=FFTExp(q+1);
    return;
end

if abs(r/Lambda)>(1-FloatingError)
    val=FFTExp(q+2);
    return;
end

pos=q+1;
if r> (Lambda/2)
    pos=pos+1;
end
del=x-D(pos);
val=FFTExp(pos);
c=1;
for i=1:M-1
    c=c*del/i;
    val=val+c*getDerAtGrid(pos,i);
end

function val=getDerAtGrid(pos,n)
global D FFTF
val=0;
for j=1:n
    val=val+FFTF(j, pos)* nchoosek(n-1,j-1);
end
val=val*exp(D(pos));
--------------------------------------
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
---------------------------------------------------------
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

--------------------------

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
%Aones=ones(1,N);
%Aones(2:2:end)=-ones(1,N/2);
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
----------------------------------------------

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
-----------------------------------------------
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
--------------------------------------
function val=GetNormPolyInt(up,low,deg) %\int^up_low x^deg e^{-x^2/2}/(2pi)^0.5 dx

if deg==0
    val=normcdf(up)-normcdf(low);
elseif deg==1
    val=-(exp(-0.5*up.^2)-exp(-0.5*low.^2))/(2*pi)^0.5;
else
    val=-(up.^(deg-1).*exp(-0.5*up.^2)-low.^(deg-1).*exp(-0.5*low.^2))/(2*pi)^0.5+...
        (deg-1)*GetNormPolyInt(up, low, deg-2);
end
--------------------------------------
%   val=\int^{D(pos+1)}_{y+D(pos)}(x-D(pos))^(deg)q(x)dx
%   deg=0,1,2,3,...,Deg, pos=1,2,...,N/2
%   
function val=GetPartInt(y, pos, deg) 
global N Deg D Left Right Lambda FFTF  PInt FloatingError
if y<0 | y >Lambda
    error('wrong y value in GetPartInt, %f', y);
end
val=PInt(deg+1, pos);
if abs(y)<FloatingError
    return;
end
d1=0;
d2=0;
d3=0;
d4=0;
d5=0;
if deg==0
    d1=-GetDerivative(D(pos), 1);
    d2=-GetDerivative(D(pos), 2);
    d3=-GetDerivative(D(pos), 3);
    d4=-GetDerivative(D(pos), 4);
    d5=-GetDerivative(D(pos), 5);
    d6=-GetDerivative(D(pos), 6);
    d7=-GetDerivative(D(pos), 7);
    val=val+d1*y+d2*y^2/2+d3*y^3/6+d4*y^4/24+d5*y^5/120+d6*y^6/720+d7*y^7/5040;
elseif deg==1
    d1=0;
    d2=-GetDerivative(D(pos),1);
    d3=-2*GetDerivative(D(pos),2);
    d4=-3*GetDerivative(D(pos),3);
    d5=-4*GetDerivative(D(pos),4);
    d6=-5*GetDerivative(D(pos),5);
    d7=-6*GetDerivative(D(pos),6);
    d8=-7*GetDerivative(D(pos),7);
    val=val+d1*y+d2*y^2/2+d3*y^3/6+d4*y^4/24+d5*y^5/120+d6*y^6/720+d7*y^7/5040+d8*y^8/40320;
elseif deg==2
    d1=0;
    d2=0;
    d3=-2*GetDerivative(D(pos),1);
    d4=-6*GetDerivative(D(pos),2);
    d5=-12*GetDerivative(D(pos),3);
    d6=-20*GetDerivative(D(pos),4);
    d7=-30*GetDerivative(D(pos),5);
    d8=-42*GetDerivative(D(pos),6);
    d9=-56*GetDerivative(D(pos),7);
    val=val+d1*y+d2*y^2/2+d3*y^3/6+d4*y^4/24+d5*y^5/120+d6*y^6/720+d7*y^7/5040+d8*y^8/40320+d9*y^9/362880;
else
    error('in GetPartInt, not implemented for deg is big than 2');
end
----------------------------------------
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
----------------------------

function void =testDensity();
global Sigma Rate Delta D N  M  Left Right  %for testing purpose
mu=(Rate-0.5*Sigma^2)*Delta;
std=Delta^0.5*Sigma;

val0=exp(-0.5*((D-mu)/std).^2)/((2*pi)^0.5*std);
val1=-val0.*(D-mu)/std^2;
val2=-val1.*(D-mu)/std^2-val0/std^2;
val3=-val2.*(D-mu)/std^2-2*val1/std^2;
y=zeros(4,N/2+1);
for i=1:(N/2+1)
y(1,i)=GetDensity(D(i),1);
y(2,i)=GetDensity(D(i),2);
y(3,i)=GetDensity(D(i),3);
y(4,i)=GetDensity(D(i),4);
end

[max(abs(y(1,:)-val0)),max(abs(y(2,:)-val1)),max(abs(y(3,:)-val2)),max(abs(y(4,:)-val3))]
exit('testDensity')

--------------------------------

function void =testDerivative();
global Sigma Rate Delta D N Left Right  %for testing purpose
mu=(Rate-0.5*Sigma^2)*Delta;
std=Delta^0.5*Sigma;
%x=D(2:N/2);%+D(1:N/2-1))/2;
x=(D(2:N/2)+D(1:N/2-1))/2;
xx=(x-mu)/std;

%((D(2:N/2)-mu)/std+(D(1:N/2-1)-mu)/std)/2;

val0=normcdf(xx);
y0=zeros(1,N/2-1);
for i=1:N/2-1;
    y0(i)=GetDerivative(x(i), 0);
end
max(abs(y0-val0))


val1= exp(-0.5*((x-mu)/std).^2)/((2*pi)^0.5*std);
y1=zeros(1,N/2-1);
for i=1:N/2-1;
    y1(i)=GetDerivative(x(i), 1);
end
max(abs(y1-val1))

val2=-val1.*(x-mu)/std^2;
y2=zeros(1,N/2-1);
for i=1:N/2-1;
    y2(i)=GetDerivative(x(i), 2);
end
max(abs(y2-val2))

val3=-val2.*(x-mu)/std^2-val1/std^2;
y3=zeros(1,N/2-1);
for i=1:N/2-1;
    y3(i)=GetDerivative(x(i), 3);
end
max(abs(y3-val3))

val4=-val3.*(x-mu)/std^2-2*val2/std^2;
y4=zeros(1,N/2-1);
for i=1:N/2-1;
    y4(i)=GetDerivative(x(i), 4);
end
max(abs(y4-val4))

%[-1, 1] 
%N=2^16
%1.110223024625157e-015
%3.730349362740526e-014
%2.650324404385174e-012
%4.292814992368221e-010
%3.808084875345230e-005

%N=2^15
%1.110223024625157e-015
%3.730349362740526e-014
%2.643218977027573e-012
%7.610651664435864e-009
%6.091333925724030e-004
%N=2^14
%1.110223024625157e-015
%3.552713678800501e-014
%6.707523425575346e-012
%2.381129888817668e-007
%0.00974608212709
%N=2^13
%1.110223024625157e-015
%3.552713678800501e-014
%3.116156221949495e-010
%7.613452908117324e-006
%0.15591597370803

%N=2^12
%9.992007221626409e-016
%1.360689338980592e-012
%1.982584763027262e-008
%2.435865389998071e-004
%2.49413345474750

%N=2^11
%2.242650509742816e-014
%1.769357993453014e-010
%1.268410187549307e-006
%0.00779339084693
%39.90301892906427

%N=2^10
%5.530242930262830e-012
%2.264817311470324e-008
%8.116278058878379e-005
%0.24928313524288
%6.379222038090229e+002

%N=2^9
%1.413602485067500e-009
%2.893335754095006e-006
%0.00518017831280
%7.94535282766810
%1.014546671040356e+004
----------------------------
%int^x_Left e^(h*u) q(u)du   for x\in [Left Right]

%test 
function void = testGetExpF(h)
global Sigma Rate Delta D N  M  Left Right  %for testing purpose
%mu=(Rate-0.5*Sigma^2)*Delta
%std=Delta^0.5*Sigma
xvalue=(D(1:(N/2-1))+D(2:(N/2)))/2;
%xvalue=D(1:(N/2-1));
yvalue=zeros(1,N/2-1);
for i=1:(N/2-1)
    yvalue(i)=GetExpF(xvalue(i), h);
end
v=exp(h*mu+0.5*std^2*h^2)*(normcdf((xvalue-mu)/std-h*std)-normcdf((Left-mu)/std-h*std));
max(abs(v-yvalue))

---------------------------

function void =testPartInt(delta);
global Sigma Rate Delta D N Left Right Lambda PInt%for testing purpose
mu=(Rate-0.5*Sigma^2)*Delta;
std=Delta^0.5*Sigma;
delta=delta*Lambda;
up=D(2:N/2);
d=D(1:N/2-1);%up;%-Lambda/2;
low=d+delta;
x=zeros(1,N/2-1);
y=zeros(1,N/2-1);

for deg=1:3
    for i=1:N/2-1
        if deg==1
            x(i)=GetNormPolyInt((up(i)-mu)/std,(low(i)-mu)/std,0);
            %x(i)=GetDerivative(up(i), 0)-GetDerivative(low(i), 0);
        elseif deg==2
            x(i)=std*GetNormPolyInt((up(i)-mu)/std,(low(i)-mu)/std,1)+...
                (mu-d(i))*GetNormPolyInt((up(i)-mu)/std,(low(i)-mu)/std,0);
        elseif deg==3
            x(i)=std^2*GetNormPolyInt((up(i)-mu)/std,(low(i)-mu)/std,2)+...
                2*std*(mu-d(i))*GetNormPolyInt((up(i)-mu)/std,(low(i)-mu)/std,1)+...
                (mu-d(i))^2*GetNormPolyInt((up(i)-mu)/std,(low(i)-mu)/std,0);
        end
        y(i)=GetPartInt(delta, i, deg-1);
    end
    [max(abs(x)),max(abs(y))]
    max(abs(x-y))
end

---------------------function void =Update(Index)  %timeIndex=NumPeriods-1:-1:0
global NumPeriods Euro Coef NumFun EpsilonContinuation N SStar KStar CP D Deg Skip Lambda Boundary
global Points Values Index  
global bug listSStar B

info=zeros(1,4);
info(1)=B;
info(2:4)=Continuation(D(info(1)));
bd=info;
flag=1;
tempB=B;

%find contination values on grid
tic;
if info(2)<=GetPayOff(D(info(1)));
    flag=0; %use same grid bounary  for previous period
end
while(flag==1)
    %Index
    %disp('temp update')
    info(1)=info-1;
    info(2:4)=Continuation(D(info(1)));
    bd=[info;bd];
    if info(2)<=GetPayOff(D(info(1)))
        tempB=info(1);
        flag=0;
    end 
end
flag=1;
info=bd(end,:);
while(flag==1)
    info(1)=info+1;
    info(2:4)=ContinuationNew(D(info(1)));
    bd=[bd;info];
    if info(2)<EpsilonContinuation;
        flag=0;
    end
    %[info(1), N/2+1]
    if info(1)>=(N/2+1)
        flag=0;
        disp('WARNING: EpsilonContinuation might be two small, in UpdateNew');
    end
end
%t=toc;
%disp('time spent for contiuation is ');
%[Index, t]
%tic;
[rows,cols]=size(bd);
coef=zeros(rows,3);
coef(:,3)=bd(:,2);
coef(:,2)=bd(:,3);
coef(:,1)=bd(:,4)/2;
disp('checking in update, the following vector should close to zeros');
coef(rows,:)
coef(rows,:)=0; 

%tic;
%SStar=FindSStar(coef(1,:),bd(1,:),bd(2,:));
KStar=Newton(bd(1,:),bd(2,:));
%t=toc;
%disp('time spent for contiuation searching sstar');
%t

Coef=coef; %update here, we might need old Coef to find KStar;
listSStar(Index+1)=SStar;
KStar=log(SStar);
[rows, void]=size(bd);
NumFun=floor((rows-2)/Deg)+1;
ind=1:Deg:(1+Deg*(NumFun-1));
col1=bd(2:end,1)';
Boundary=col1(ind);
Boundary=[bd(1,1),Boundary];
%disp('leaving update')
%[Index,SStar, KStar, NumFun, Index]
%Coef
%Boundary
%exit(1)

if Index==0
    listSStar'
end
--------------------------

%function void = main(model, strike, fftpower,  skip1,skip2, deg,m, left, right)

function void = main(model, strike, fftpower, deg,m, left, right)

global KStar NumPeriods ScaleUnit Strike Lambda  Coef Deg Boundary NumFun Index
global bug M D% for debug, remove late 
%global Skip1 Skip2
global Skip
%T=0
%for i=1:10000
%tic;

%Init(model, strike, fftpower, skip1, skip2,  deg,m, left, right);
Constant(model);
Init(model, strike, fftpower, deg,m, left, right);

%European();

%test 
%testPInt(model);
%exit('tempin main');

%test
testPartInt(strike)
exit('temp exit in main')

%test
%testDerivative()
%exit('temp exit in main')

for Index=NumPeriods-1:-1:0
    bug=NumPeriods-2;
    UpdateNew();
end
%exit('temp exit in main');
%x=[36,38,40,42,44]';
x=100;
[rows, cols]=size(x);
y=zeros(rows,1);
c=zeros(rows,1);
logStock=log(x/ScaleUnit);
exp(logStock)*ScaleUnit
counter=1;
for k=1:rows
    if logStock(k)<KStar
        c(k)=GetPayOff(logStock(k));
    else
        for m=1:NumFun
            if D(Boundary(m))>=logStock(k)
                c(k)=polyval(Coef(m,:),exp(logStock(k)));
                break;
            else
                ;
            end
        end
        if m==NumFun
            disp('not find the value for ')
            exp(logStock(k))*ScaleUnit   
        end            
    end
end
ScaleUnit*c
%t2=toc;
%T=T+t2;
%end
%T/10000=0.17623696944723

--------------------------
%   PInt(deg, j)=\int^{D(j+1)}_{D(j)}(x-D(j))^(deg-1)q(x)dx
%   deg=1,2,3,...,Deg+1, j=1,2,...,N/2
%   \int^{D(N/2+1)}_{D(N/2))(x-D(j))^(deg-1)q(x)dx is set to be 0 for all deg=1,2,Deg+1

function InitFFT(model) 
global N Deg D Left Right Lambda FFTF FFTExp CDF PInt 

bma=Right-Left;
bmapi=bma/pi;
Aones=ones(1,N);
Aones(2:2:end)=-ones(1,N/2);
K=0:1:N-1;
%Y=0:1:(N-1);
%Y=2*bma*Y/N+2*Left-Right; %same as Y*Lambda+2*Left-Right

F=DensityCoeficient(model, Left,Right, N); 

FFTF=zeros(M,N); % Derivative of q(x) up to order M-1, M has to be odd!
F1=Aones.*F;
FFTF(1,:)=N*real(ifft(F1));  % order 0 derivative q value
for i=1:((M-1)/2)
    F1=F1.*K*pi/(Right-Left);  % 
    FFTF(2*i,:)=(-1)^i*N*imag(ifft(F1)); % order 2*i-1 derivative
    F1=F1.*K*pi/(Right-Left);  %  
    FFTF(2*i+1,:)=(-1)^i*N*real(ifft(F1)); % order 2*i derivative
end
FFTF=FFTF(:,N/2+1:N);

%init FFTAD
FFTAD=zeros(Deg+1,N);
F1=F(2:N)./K(2:N); %first order anti, i.e. CDF
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

FFTAD=FFTAD(:,N/2+1:N);
CDF=FFTAD(1,:);

PInt=zeros(Deg+1,N/2-1);
for deg=1:Deg+1
    PInt(deg,:)=FFTAD(1,2:N/2)*Lambda^(deg-1)-FFTAD(1,1:N/2-1)*0^(deg-1);
    factor=1;
    for k=1:deg-1
        factor=factor*(deg-k) ;
        PInt(deg,:)=PInt(deg,:)+(-1)^k*factor*(FFTAD(k+1,2:N/2)*Lambda^(deg-k-1)-...
            FFTAD(k+1,1:(N/2-1))*0^(deg-k-1) );
    end
end
PInt=[PInt, zeros(Deg+1,1)];

---------------------------
%function void = Init(modelType, strike, fftpower, skip1,skip2, deg,m, left, right) 

function void = Constant(modelType) 

%constant 
global EpsilonSearching EpsilonContinuation Eps  Discount
global Maturity Sigma Rate NumPeriods Delta  Oalpha Obeta
%model paramters 
global Mq Mlambda Mmu Mgamma Mdelta %Merton's Model
global Kp Klam Klamp Klamm %Kou's model
global NIGmu NIGdelta NIGalpha NIGbeta %Normal Inverse Gaussian Model
global VGnu VGmu VGtheta %VG model
global TSgamma TSalpha TSalphap TSalpham TSlambdap TSlambdam TScp TScm %Tempered Stable Process

Oalpha=1;
Obeta=1; 
EpsilonSearching=10^(-15);
EpsilonContinuation=10^(-15);
FloatingError=10^(-15);
%Rate=0.5;
Rate=0.1;
%Rate=0.06;
%NumPeriods=3;
NumPeriods=50;
Eps=(0:1:NumPeriods)*(10^(-4)/NumPeriods);
Maturity=1;  
Delta=Maturity/NumPeriods;
Discount=exp(-Rate*Delta);

------------------------------------------------
    

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
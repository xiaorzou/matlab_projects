%x do not have to be in grid
%S=exp(x): the price of stock
%val=continuation value
%d1: first derivative
%d2: second derivative
function vector=ContinuationNew(x)
global KStar Oalpha Obeta Strike Boundary Discount NumFun Coef NumPeriods EpsilonSearching Lambda
global CP D Deg Euro
global Index Points KP Left Right  SStar %donot need this one late.

if x<Left | x>Right
    error('wrong range for x in contiaiton');
end
val=0;
s=exp(x);

low=min(max(KStar-x, Left), Right);
up1=min(max(log(Strike)-x, Left), Right); %log(Strike)=0


q1low=GetDensity(low,1); %check
%q1up1=GetDensity(up1,1);
q2low=GetDensity(low,2); %check

e0low=GetExpF(low,0); %check
e0up1=GetExpF(up1,0); %check

e1low=GetExpF(low,1); %check
e1up1=GetExpF(up1,1); %check

up2=min(max(D(Boundary(2))-x, Left), Right);
part2=Oalpha*Strike*(e0up1-e0low)-Obeta*s*(e1up1-e1low);
d1=q1low*(Obeta*SStar-Oalpha*Strike)/s-Obeta*e1low+polyval(Coef(1,:), SStar)*q1low/s;
d2=Oalpha*Strike*q1low/s^2-(Oalpha*Strike-Obeta*SStar)*q2low/s^2;

d2=d2-polyval(Coef(1,:), SStar)*(q1low+q2low)/s^2+polyval((Coef(1,:).*(Deg:-1:0)), SStar)*q1low/s^2;

%get part3, compute int^(infity)_(KStar-x)C(exp(u+x))q(u)du
part3=0;

%First, deal with\int^{D(Boundary(2))-x}_{KStar-x}C(exp(u+x))q(u)du
cc=0; 
for i=1:(Deg+1)
    temp=s^(Deg+1-i)*( GetExpF(up2,Deg+1-i)-GetExpF(low,Deg+1-i));
    part3=part3+Coef(1,i)*temp;
    cc=cc+(Deg+1-i)*Coef(1,i)*temp; 
    d2=d2+(Deg+1-i)^2*Coef(1,i)*temp; 
end

if Index==100
    [3, cc, d2]
end
if part3<0 
   % part3
end

%deal with \int^{D(Boundary(j+1))}_{D(Boundary(j))}
for j=2:NumFun
    lowt=min(max(D(Boundary(j))-x, Left), Right);
    upt=min(max(D(Boundary(j+1))-x, Left), Right);
    for i=1:(Deg+1)
      temp=s^(Deg+1-i)*( GetExpF(upt,Deg+1-i)-GetExpF(lowt,Deg+1-i) );
      part3=part3+Coef(j,i)*temp;
      cc=cc+(Deg+1-i)*Coef(j,i)*temp;
      d2=d2+(Deg+1-i)^2*Coef(j,i)*temp;      
      if Index==100
          [i, s,GetExpF(upt,Deg+1-i), GetExpF(lowt,Deg+1-i), Coef(i,i)]
      end
    end 
end

d1=Discount*(d1+cc/s);
d2=Discount*(d2-cc/s^2);
%val=Discount*(-part2+part3+Oalpha*Strike*e0up2-Obeta*s*e1up2);
val=Discount*(-part2+part3+Oalpha*Strike*e0up1-Obeta*s*e1up1);
%disp('leave continuation')
vector=[max(val,0), d1, d2];
if val <0
    val
end

************************
function coefdensity=DensityCoeficient(modelType, leftend,rightend, N)
x=0:1:(N-1);
coefdensity=MyChar(modelType, x*pi/(rightend-leftend)).*exp(-x*leftend*pi/(rightend-leftend)*complex(0,1));
coefdensity=2*real (coefdensity)/(rightend-leftend);
coefdensity(1)=coefdensity(1)/2;

******************************
%q_{(order)}(at): if order=0, the value of q(at), otherwise the value of
%antiderivative, expend around D(lpos)
function val=EvaluateDensity(at, lpos, order)
global D FFT FFT3 Lambda
delta=at-D(lpos);
if abs(delta)>Lambda
    error('wrong in EvaluateDensity')
end
if order==0
    if delta==0
        val=FFT(3, lpos);
    else
        val=FFT(3, lpos)+delta*FFT(2,lpos)+0.5*delta^2*FFT(1,lpos)+delta^3*FFT3(1,lpos)/6+delta^4*FFT3(2, lpos)/24 +delta^5*FFT3(3, lpos)/120 +delta^6*FFT3(4, lpos)/720;
    end
elseif order ==1 % \int^at_{-\infty} q(x)dx
    if delta==0
        val=FFT(4,lpos);
    else
        val=FFT(4,lpos)+delta*FFT(3,lpos)+delta^2*FFT(2,lpos)/2+delta^3*FFT(1,lpos)/6+delta^4*FFT3(1, lpos)/24 +...
        delta^5*FFT3(2, lpos)/120+delta^6*FFT3(3, lpos)/720+delta^7*FFT3(4, lpos)/5040;
    end
elseif order==2
    if delta==0
        val=FFT(5,lpos);
    else
        val=FFT(5,lpos)+delta*FFT(4,lpos)+delta^2*FFT(3,lpos)/2+delta^3*FFT(2,lpos)/6 + ...
        delta^4*FFT(1, lpos)/24 +delta^5*FFT3(1, lpos)/120 +...
        delta^6*FFT3(2, lpos)/720+delta^7*FFT3(3, lpos)/5040+delta^8*FFT3(4, lpos)/40320;
    end
elseif order==3
    if delta==0
        val=FFT(6,lpos);
    else
        val=FFT(6,lpos)+delta*FFT(5,lpos)+delta^2*FFT(4,lpos)/2+delta^3*FFT(3,lpos)/6 + ...
        delta^4*FFT(2, lpos)/24 +delta^5*FFT(1, lpos)/120 +...
        delta^6*FFT3(1, lpos)/720+delta^7*FFT3(2, lpos)/5040+delta^8*FFT3(3, lpos)/40320+...
        delta^9*FFT3(4, lpos)/362880;
    end
elseif order==4
    if delta==0
        val=FFT(7,lpos);
    else
        val=FFT(7,lpos)+delta*FFT(6,lpos)+delta^2*FFT(5,lpos)/2+delta^3*FFT(4,lpos)/6 + ...
        delta^4*FFT(3, lpos)/24 +delta^5*FFT(2, lpos)/120 +...
        delta^6*FFT(1, lpos)/720+delta^7*FFT3(1, lpos)/5040+delta^8*FFT3(2, lpos)/40320+...
        delta^9*FFT3(3, lpos)/362880+delta^10*FFT3(4, lpos)/3628800;
    end
elseif order==5
    if delta==0
        val=FFT(8,lpos);
    else
        val=FFT(8,lpos)+delta*FFT(7,lpos)+delta^2*FFT(6,lpos)/2+delta^3*FFT(5,lpos)/6 + ...
        delta^4*FFT(4, lpos)/24 +delta^5*FFT(3, lpos)/120 +...
        delta^6*FFT(2, lpos)/720+delta^7*FFT(1, lpos)/5040+delta^8*FFT3(1, lpos)/40320+...
        delta^9*FFT3(2, lpos)/362880+ delta^10*FFT3(3, lpos)/3628800 +delta^11*FFT3(4, lpos)/39916800;
    end
else
    error('not implemented yet, in EvaluateDensity');
end
*****************************************************

function coef=FindCoef(bd)
global Deg D Index Coe

if Deg~=4
    error('Deg has to be 4!, in FindCoef');
end
[rows, void]=size(bd);
numfun=floor((rows-2)/Deg)+1;
coef=zeros(numfun,Deg+1);

s=[exp(D(bd(1,1))),exp(D(bd(2,1)))];
p=s(2)-s(1);

Coe(5)=bd(1,2);
Coe(4)=bd(1,3)*p;
Coe(3)=bd(1,4)*p^2/2;
Coe(1)=(bd(2,4)/6+bd(1,4)/3)*p^2+bd(1,3)*p+bd(1,2)-bd(2,2);
Coe(2)=bd(2,2)-p^2*bd(1,4)/2-bd(1,3)*p-bd(1,2)-Coe(1);


coef(1,1)=Coe(1)/p^4;
coef(1,2)=-4*Coe(1)*s(1)/p^4+Coe(2)/p^3;
coef(1,3)=6*s(1)^2*Coe(1)/p^4-3*s(1)*Coe(2)/p^3+Coe(3)/p^2;
coef(1,4)=-4*s(1)^3*Coe(1)/p^4+3*s(1)^2*Coe(2)/p^3-2*s(1)*Coe(3)/p^2+Coe(4)/p;
coef(1,5)=Coe(1)*(s(1)/p)^4-Coe(2)*(s(1)/p)^3+Coe(3)*(s(1)/p)^2-Coe(4)*s(1)/p+Coe(5);

e=bd(1,2);
d=bd(1,3);
c=bd(1,4)/2;

a=((bd(2,4)-2*c)/p-6*(bd(2,2)-c*p^2-d*p-e)/p^3)/(6*p);
b=(bd(2,2)-c*p^2-d*p-e)/p^3-p*a;
coef=zeros(1,5);
coef(1,1)=a;
coef(1,2)=-4*a*s(1)+b;
coef(1,3)=6*s(1)^2*a-3*s(1)*b+c;
coef(1,4)=-4*s(1)^3*a+3*s(1)^2*b-2*s(1)*c+d;
coef(1,5)=a*s(1)^4-b*s(1)^3+c*s(1)^2-d*s(1)+e;

for i=1:(numfun-1)
    s=exp(D(bd( (2+(i-1)*Deg):(2+i*Deg),1)))';
    y=bd((2+(i-1)*Deg):(2+i*Deg),2);
    coef(i+1,:)=polyfit(s, y, Deg);
end

************************************
function sstar=FindSStar(coef,info1,info2)
global Strike Oalpha Obeta D Deg Index Coe
payoff1=GetPayOff(D(info1(1)));
payoff2=GetPayOff(D(info2(1)));
%[payoff1, payoff2]
%[info1(2), info2(2)]
s=[exp(D(info1(1))), exp(D(info2(1)))];
sstar=0;
if payoff1<info1(2) | payoff2>info2(2)
    [payoff1-info1(2), payoff2-info2(2)]
    error('wrong in FindSStar, can not find SStar');
end
coef(4)=coef(4)+Obeta;
coef(5)=coef(5)-Oalpha*Strike;
sstar=NewtonNew(coef, Deg, s(1), s(2));
*******************************************
function sstar=FindSStarNew(info1,info2)
global Strike Oalpha Obeta D Deg Index Coe 
payoff1=GetPayOff(D(info1(1)));
payoff2=GetPayOff(D(info2(1)));
s=[exp(D(info1(1))), exp(D(info2(1)))];
sstar=0;

coef=Coe;
coef(4)=coef(4)+Obeta*(s(2)-s(1));
coef(5)=coef(5)+Obeta*s(1)-Oalpha*Strike;

sstar=NewtonNew(coef, Deg, 0, 1);
if sstar<0 | sstar>1
    error('wrong in FindSStarNew')
end
sstar=sstar*(s(2)-s(1))+s(1);
*****************************************
%q(x)^(order)
function val=GetDensity(x, order) %order-1 is the degree of derivative, order=1, ..., M-1
global D FFTF Lambda M Left Right N

if order>=M
    error('in GetDensity, the order %d is bigger than M-1:%d ', order, M-1);
end
if (x<Left | x>Right) %& order==0 , for order >0, if x is out of range, what we should return?
    val=0;
    return;
end
q=floor((x-Left)/Lambda);
r=x-Left - Lambda *q;
if q>=N/2
    q=N/2-1;
    r=0;
end
if (r/Lambda)<10^(-14)
%    disp('enter Gend1');
    val=FFTF(order,q+1);
    %[order,val]
    return;
elseif (r/Lambda)>(1-10^(-14))
    val=FFTF(order, q+2);
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
    val=FFTF(order,pos);
    for i=1:M-order
        c=c*delta/i;
        val=val+c*FFTF(i+order,pos);
    end
end
***************************************************
%int^x_Left e^(h*u) q(u)du   for x\in [Left Right]
function val = GetExpF(x, h)
global N Left Right Lambda  FFTF FFTG D M FFTAD

if x<Left | x>Right
    error('Wrong range for x %d', x)
end

q=floor((x-Left)/Lambda);
r=x-Left - Lambda *q;

if q==N/2
    q=q-1;
end

if abs(r/Lambda)<10^(-12)
    if h==0
        val=FFTAD(1, q+1);
    else
        val=FFTG(h, q+1);
    end
    return;
end

if abs(r/Lambda)>(1-10^(-12))
    if h==0
        val=FFTAD(1, q+2);
    else
        val=FFTG(h, q+2);
    end
    return;
end

pos=q+1;
if r> (Lambda/2)
    pos=pos+1;
end

del=x-D(pos);

if h==0
    val=FFTAD(1,pos);
else
    val=FFTG(h,pos);
end
%val=val+del*exp(h*D(pos))*FFTF(1,pos)+del^2*exp(h*D(pos))*( h*FFTF(1,pos)+FFTF(2,pos))/2+...
%    del^3*exp(h*D(pos))*( h^2*FFTF(1,pos)+2*h*FFTF(2,pos) + FFTF(3,pos))/6+...
%    del^4*exp(h*D(pos))*( h^3*FFTF(1,pos)+3*h^2*FFTF(2,pos) + 3*h*FFTF(3,pos)+ FFTF(4,pos))/24+...
%    del^5*exp(h*D(pos))*( h^4*FFTF(1,pos)+4*h^3*FFTF(2,pos) + 6*h^2*FFTF(3,pos)+ 4*h*FFTF(4,pos) + FFTF(5,pos))/120+...
%    del^6*exp(h*D(pos))*( h^5*FFTF(1,pos)+5*h^4*FFTF(2,pos) +
%    10*h^3*FFTF(3,pos)+ 10*h^2*FFTF(4,pos) + 5*h*FFTF(5,pos)+FFTF(6,pos))/720;
c=1;
for i=1:M-1
    c=c*del/i;
    val=val+c*getDerAtGrid(pos,i,h);
end

function val=getDerAtGrid(pos,n,h)
global D FFTF
val=0;
for j=1:n
    val=val+h^(n-j)*FFTF(j, pos)* nchoosek(n-1,j-1);
end
val=val*exp(h*D(pos));
***********************************
%F(0,y)=q(y)
%F(1,y)=\int^{y}_{Left}q(x)dx,
%F(2,y) meet the condtion F(s,y)'=F(s-1,y) and F(s,0)=0
%s=0, 1, 2, 3, 4

function val=GetF(s,y) 


global  N F Left Right
%if y<Left
%    y=Left;
%elseif y>Right
%    y=Right;
%end

if y<Left | y>Right
    disp('y is not in the range [Low, Up]')
    [Left, y, Right]
    error('wrong range for y value in GetF')
end

piobma=pi/(Right-Left);
kk=0:1:(N-1);
val=0;
if s==0
    val=sum(F*cos(kk*piobma*(y-Left)));
elseif s==1
    val=F(1)*(y-Left)+sum(F(2:N).*sin(kk(2:N)*piobma*(y-Left))./kk(2:N))/piobma;
elseif s==2
    val=0.5*F(1)*(y-Left)^2+ sum(F(2:N).*(1-cos(kk(2:N)*piobma*(y-Left)))./(kk(2:N).^2))/piobma^2;
elseif s==3
    val=(y-Left)/piobma^2*sum(F(2:N)./(kk(2:N).^2))-sum(sin(kk(2:N)*piobma*(y-Left)).*F(2:N)./(kk(2:N).^3))/piobma^3+F(1)*(y-Left)^3/6;
elseif s==4
    val=F(1)*(y-Left)^4/24+0.5*(y-Left)^2*sum(F(2:N)./(kk(2:N).^2))/piobma^2-sum(F(2:N).*(1- cos(kk(2:N)*piobma*(y-Left)))./(kk(2:N).^4))/piobma^4;
elseif s==5
    val=-(y-Left)/piobma^4*sum(F(2:N)./(kk(2:N).^4))+(y-Left)^3/piobma^2*sum(F(2:N)./(kk(2:N).^2))/6+sum(sin(kk(2:N)*piobma*(y-Left)).*F(2:N)./(kk(2:N).^5))/piobma^5+F(1)*(y-Left)^5/120;
else
    error('not yet implmemt! in GetF');
end
****************************************************
function void =GetFFT(F)
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

%mu=(Rate-0.5*Sigma^2)*Delta
%std=Delta^0.5*Sigma
%(Left-mu)/std-h*std
%for h=1:Deg
%    v=exp(h*mu+0.5*std^2*h^2)*(normcdf((D(1:N/2)-mu)/std-h*std)-normcdf((Left-mu)/std-h*std));
%    max(v-FFTG(h,:))
%end

%exit('tme exit');

%%%%%%%%%%%%%%%%%  do not remove below, can be used in future    %%%%%%%%%%%%%%%%%%%%%%%

%FFT for antiderivative, we do not need this right now.
FFTAD=zeros(1,N);
F1=F(2:N)./K(2:N);
F1=[0,F1];
F1=F1.*Aones;
FFTAD(1,:)=F(1)*(Y-Left)+N*bma*imag(ifft(F1))/pi;
FFTAD=FFTAD(N/2+1:N);

%second order anti
%F1(2:N)=F1(2:N)./K(2:N);
%temp2=N*real(ifft(F1));
%FFTAD(2,:)=0.5*F(1)*(Y-Left).^2-bmapi^2.*temp2+ temp2(N/2+1)*bmapi^2;

%third order
%F1(2:N)=F1(2:N)./K(2:N);
%FFTAD(3,:)=(1/6)*F(1)*(Y-Left).^3-bmapi^3*(N*imag(ifft(F1)))+(Y-Left)*bmapi^2*temp2(N/2+1);

%fourth order
%F1(2:N)=F1(2:N)./K(2:N);
%temp3=N*real(ifft(F1));
%FFTAD(4,:)=(F(1)/24)*(Y-Left).^4+bmapi^4*temp3-bmapi^4*temp3(N/2+1)+0.5*(Y-Left).^2*bmapi^2*temp2(N/2+1);

%fifth order
%F1(2:N)=F1(2:N)./K(2:N);
%temp4=N*imag(ifft(F1));
%FFTAD(5,:)=(F(1)/120)*(Y-Left).^5+bmapi^5*temp4-bmapi^4*temp3(N/2+1)*(Y-Left)+bmapi^2*temp2(N/2+1)*(Y-Left).^3/6;

%FFTAD=FFTAD(:,N/2+1:N);

******************************************
%F(0,y)=q(y)
%F(1,y)=\int^{y}_{Left}q(x)dx,
%F(2,y) meet the condtion F(s,y)'=F(s-1,y) and F(s,0)=0
%s=0, 1, 2, 3, 4
function val=GetFNew(s,y) 
global  N Left Right Lambda

if y<Left | y>Right
    disp('y is not in the range [Left, Right]')
    [Left, y, Right]
    error('wrong range for y value in GetFNew')
end

lpos=floor((y-Left)/Lambda)+1;
lpos=min(max(lpos,1),N/2);
val=EvaluateDensity(y, lpos, s);

********************************************
function val=GetPartInt(coef,low,up,y) %\int^{up-y}_{low-y}coef(x+y)*q(x)dx
global  Deg tempIndex NumPeriods F Right Left Boundary KStar Points KP

bug=1000;
if low>up %to be removed for efficiency
    [low up, KP]
    [Boundary(1), KStar, Boundary(2), Boundary(3)]
    Points
    error('wrong range for low and up in GetPartInt');
end

if up-y>Right
    up=Right+y;
end
if low-y<Left
    low=y+Left;
end

val=polyval(coef,up)*GetF(1,up-y)-polyval(coef,low)*GetF(1,low-y);
coef=coef(1:Deg).*(Deg:-1:1);
val=val-(polyval(coef,up)*GetF(2,up-y)-polyval(coef,low)*GetF(2,low-y));
if Deg>=2
    coef=coef(1:(Deg-1)).*((Deg-1):-1:1);
    val=val+(polyval(coef,up)*GetF(3,up-y)-polyval(coef,low)*GetF(3,low-y));
    if Deg>=3
        coef=coef(1:(Deg-2)).*((Deg-2):-1:1);
        val=val-(polyval(coef,up)*GetF(4,up-y)-polyval(coef,low)*GetF(4,low-y));
        if Deg>=4
            coef=coef(1:(Deg-3)).*((Deg-3):-1:1);
            val=val+(polyval(coef,up)*GetF(5, up-y)-polyval(coef,low)*GetF(5,low-y));
            if Deg>=5
                error('not implement yet, n GetPartInt');
            end
        end
    end
end
if tempIndex==NumPeriods-2+bug
val
exit('leave getpartint')
end

***********************************************

function val=GetPartIntNew(coef,low,up,y) %\int^{up}_{low}coef(x+y)*q(x)dx
global  Deg tempIndex NumPeriods F Right Left Boundary KStar Points KP

bug=1000;
if low>up %to be removed for efficiency
    [low, up, KP]
    [Boundary(1), KStar, Boundary(2), Boundary(3)]
    Points
    error('wrong range for low and up in GetPartInt');
end

%if up>Right
%    up=Right;
%end
%if low<Left
%    low=Left;
%end
%[low,up]
low=min(max(low,Left),Right);
up=min(max(up,Left),Right);
%[low, up, Left, Right]
val=polyval(coef,up+y)*GetFNew(1,up)-polyval(coef,low+y)*GetFNew(1,low);
coef=coef(1:Deg).*(Deg:-1:1);
val=val-(polyval(coef,up+y)*GetFNew(2,up)-polyval(coef,low+y)*GetFNew(2,low));
if Deg>=2
    coef=coef(1:(Deg-1)).*((Deg-1):-1:1);
    val=val+(polyval(coef,up+y)*GetFNew(3,up)-polyval(coef,low+y)*GetFNew(3,low));
    if Deg>=3
        coef=coef(1:(Deg-2)).*((Deg-2):-1:1);
        val=val-(polyval(coef,up+y)*GetFNew(4,up)-polyval(coef,low+y)*GetFNew(4,low));
        if Deg>=4
            coef=coef(1:(Deg-3)).*((Deg-3):-1:1);
            val=val+(polyval(coef,up+y)*GetFNew(5, up)-polyval(coef,low+y)*GetFNew(5,low));
            if Deg>=5
                error('not implement yet, n GetPartInt');
            end
        end
    end
end
if tempIndex==bug
val
exit('leave getpartint')
end

************************************************
%function val=GetPayOff(pos)
function val=GetPayOff(logPrice)
global Oalpha Obeta Strike 
val=max(Oalpha*Strike-Obeta*exp(logPrice),0);

****************************************
%function val=MyFFT(modelType,sp,sn)
function val=MyFFT(modelType)
global N Alpha Eta Lambda  DP U FFTError
vValues=0:1:N-1;
vValues=vValues*Eta;

xValues=zeros(1,N);
%val=zeros(DP+1, sp-sn+1); % we have degree -DN,...0,...,DP
val=zeros(DP+1, N);
for j=1:DP+1
    xValues(1)=MyChar(modelType, vValues(1)-complex(0,1)*(Alpha+j-1))/(complex(0,1)*vValues(1) +Alpha)*exp(complex(0,1)*N*vValues(1)*Lambda/2) ;
    for k=2:N
        xValues(k)=MyChar(modelType, vValues(k)-complex(0,1)*(Alpha+j-1))/(complex(0,1)*vValues(k) +Alpha)*exp(complex(0,1)*N*vValues(k)*Lambda/2)*(3+(-1)^k);
    end
    temp=fft(xValues);
    %val(j,:)=(Eta/(3*pi))*(exp(-Alpha*U).*temp(sn:sp));
    val(j,:)=(Eta/(3*pi))*(exp(-Alpha*U).*temp);
end
val=max(real(val),0);
[a,ll]=min(abs(val(1,:)-1));
[b,u]=min(val(1,ll:N));
lowPos=ll;
upPos=ll+u-1;
lflag=0;
%uflag=0;
i=ll;
while (lflag==0)
    if abs(val(1,i)-1)>FFTError & lflag==0
        lowPos=i-1;
        lflag=1;
    end
    i=i+1;
end

val=val(:,lowPos:upPos);
U=U(lowPos:upPos);
N=upPos-lowPos+1;
******************************************
function void =UpdateNew(Index)  %timeIndex=NumPeriods-1:-1:0
global NumPeriods Euro Coef NumFun EpsilonContinuation N SStar KStar CP D Deg Skip Lambda Boundary
global Points Values Index  
global bug listSStar 

info=zeros(1,4);
info(1)=Boundary(1);
info(2:4)=ContinuationNew(D(info(1)));
bd=info;
flag=1;
tic;
if info(2)<=GetPayOff(D(info(1)));
    flag=0;
end
while(flag==1)
    info(1)=FindPreviousPos(info);
    info(2:4)=ContinuationNew(D(info(1)));
    bd=[info;bd];
    if info(2)<=GetPayOff(D(info(1)))
        flag=0;
    end 
end
flag=1;
info=bd(end,:);
while(flag==1)
    info(1)=FindNextPos(info);
    info(2:4)=ContinuationNew(D(info(1)));
    bd=[bd;info];
    if info(2)<EpsilonContinuation;
        flag=0;
    end
    if info(1)>=(N/2+1)
        flag=0;
        disp('WARNING: EpsilonContinuation might be two small, in UpdateNew');
    end
end
coef=FindCoef(bd);
if Index==46
    disp('check find coef')
    coef(1,:)
    [polyval(coef(1,:), exp(D(bd(1,1))))-bd(1,2), polyval(coef(1,:),exp(D(bd(2,1))))-bd(2,2)]
end
SStar=FindSStarNew(bd(1,:),bd(2,:));
Coef=coef; %update here, we might need old Coef to find KStar;
listSStar(Index+1)=SStar;
KStar=log(SStar);
[rows, void]=size(bd);
NumFun=floor((rows-2)/Deg)+1;
ind=1:Deg:(1+Deg*(NumFun-1));
col1=bd(2:end,1)';
Boundary=col1(ind);
Boundary=[bd(1,1),Boundary];
if Index==0
    listSStar'
end
%function void = Init(modelType, strike, fftpower, skip1,skip2, deg,m, left, right) 

function void = Constant(modelType, numPeriods,  alpha, beta) 

%constant 
global EpsilonSearching EpsilonContinuation Eps  Discount
global Maturity Sigma Rate NumPeriods Delta  Oalpha Obeta
%model paramters 

global ModelScaling

global Mq Mlambda Mmu Mgamma Mdelta %Merton's Model
global Kp Klam Klamp Klamm %Kou's model
global NIGmu NIGdelta NIGalpha NIGbeta %Normal Inverse Gaussian Model
global VGnu VGmu VGtheta %VG model
global TSgamma TSalpha TSalphap TSalpham TSlambdap TSlambdam TScp TScm %Tempered Stable Process


Oalpha=alpha;
%Obeta=1; 
Obeta=beta;
EpsilonSearching=10^(-15);
EpsilonContinuation=10^(-15);
FloatingError=10^(-15);
%Rate=0.5;
Rate=0.1;
%Rate=0.06;
NumPeriods=numPeriods;
%NumPeriods=50;
Eps=(0:1:NumPeriods)*(10^(-4)/NumPeriods);
Maturity=1;  
Delta=Maturity/NumPeriods;
Discount=exp(-Rate*Delta);
if modelType==1 
    %Sigma=1;
    Sigma=0.1305;          %volatility constant.
    %Sigma=0.2;
elseif modelType==2 %Merton model
    Mq=0.0;
    Sigma=0.0671654;
    Mlambda=1.61843;
    Mdelta=0.03235;
    Mmu=-0.08650;
    Mgamma=Rate-Mq-0.5*Sigma^2-Mlambda*(exp(0.5*Mdelta^2+Mmu)-1);
    %Delta=Maturity/NumPeriods;
elseif modelType==3 %Kou model
    Sigma=0.06498;
    Klam=4.13590;
    Klamp=24.22148;
    Klamm=24.22148;
    Kp=0.09006;
elseif modelType==4 %NIG
    NIGalpha=28.42141;
    NIGbeta=-15.08623;
    NIGmu=0.05851;
    NIGdelta=0.31694;
elseif modelType==5 %VG model lack other information, using cutoff=800, fftpower=12, alpha=0.5
    VGnu=0.2;
    Sigma=0.12;
    VGmu=-0.14;
    VGtheta=-0.14; %??
    %Delta=Maturity/NumPeriods;
elseif modelType==6 %TS
    TSalpha=0.273;
    TSalphap=0.273; %???
    TSalpham=0.273; %??? 
    TScp=2.093;
    TScm=1.952;
    TSlambdap=38.209;
    TSlambdam=16.050;
    %Delta=Maturity/NumPeriods; 
end
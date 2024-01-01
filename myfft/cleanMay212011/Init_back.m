%function void = Init(modelType, strike, fftpower, skip1,skip2, deg,m, left, right) 

function void = Init_back(modelType, strike, fftpower, deg,m, left, right, numA, numB) 
global ScaleUnit Lambda  NumPeriods
%data structure for the implementation 
global N NumA NumB  %nodes in FFT
%global ST
global D   %for grid points over the range [Left, Right] for W, dim(D)=N/2+1, a row vector
global Left Right %the valid range for q_W(w), q_W(w)=0 if x is outside the range. X=Scaling*W
global Deg  %approximate contination function locally using Deg-degree polynomial

%global M    %M-1: the degree of polynomial (Taylor expansion) 
            %that is used to approximate continuation at points that is not
            %on grid. M has to be odd in the inplmetation
%global FFTF %FFTF(1,pos)=q(D(pos)), FFTF(j,x)=FFTF(j-1,x)', j=2,3,..M


global FFTExp %FFTExp(i)=\int^D(i)_Left exp(Scaling*w)q_W(w)dw,

global FFTAD%FFTAD(i)=\int^D(i)_Left q_W(w)dw

%global PIntA %integration PInt(deg, j)=\int^{D(j)+DeltaA}_{D(j)}(x-D(j))^(deg-1)q(x)dx
%global PIntB %integration PInt(deg, j)=\int^{D(j)+DeltaB}_{D(j)}(x-D(j))^(deg-1)q(x)dx

global Euro %for european price
global Euro2

%global Skip1 % used to refine the searching nodes forward
%global Skip2 %used to refine the searching nodes afterward

global CP %center position of D, is equal to 1+N/4

global Fib %finnahi ? 

global D3

global Strike

%global CoefA %for coefficents of contiuation before CP:  c1(x-D(pos))^Deg+... where D(pos) is the starting point of the interval 
%global CoefB %for coefficents of contiuation after CP

global Coef

%global SP % SStar is in [BA(SP), BA(SP+1))

%global BA % the position of boundary points before CP
%global BB % the position of boundary points after CP
global B

global SA % Skip for grid points before CP
global SB % Skip for grid points after CP

%global GPA % grid points before CP
%global GPB % grid points after CP
global GP

%global NumFun  %dynamically keep the record the number of functions we have at current time.

global F Index listSStar  %temperoarly, not need in the future

global Scaling

N=2^fftpower;
CP=N/4+1; %center point position, D(CP)=0!
Scaling=1;

Left=left; %[LEFT RIGHT] is the cutoff for the standard with std=1, mu=0
Right=right; 

Deg=deg;

%M=m; 
%if floor(M/2)==M/2
%    error('Init: M has to be odd');
%end

%Left=left;
%Right=right;

%Lambda=2*(Right-Left)/N; %old
Lambda=(Right-Left)/(N/2); %new

D=(0:1:N/2)*Lambda+Left; %N/2+1 elements
Euro=zeros(1,N/2+1);

ScaleUnit=strike;        
Strike=strike/ScaleUnit; %must be 1!
%KStar=log(Strike);
%SStar=exp(KStar);

%N=2^10; %for testing 
%CP=N/4+1 %for testing 
%Deg=9; %for testing 
%Left=-1; %for testing 
%Right=1; %for testing 
%Lambda=2*(Right-Left)/N; %for testing 
%D=(0:1:N/2)*Lambda+Left %for testing 

NumA=numA; %should be adjusted to cover the exercise boundaries
%SP=NumA; %init it as the number such that exercise boundary is in [BA(NumA),BA(NumA+1)]
SA=1; %flexible 
st=CP-Deg*NumA*SA;
GP=st+(0:1:Deg*NumA)*SA;
B=st+(0:1:NumA)*SA*Deg;
%CoefA=zeros(NumA,Deg+1);

NumB=numB; %flexible 
SB=floor((N/4)/(Deg*NumB));
if SB==0
    error('error in init, need reduce the numB');
end
gbp=CP+(1:1:Deg*NumB)*SB;
GP=[GP,gbp];
bb=CP+(1:1:NumB)*Deg*SB;
B=[B,bb];

%CoefB=zeros(NumB,Deg+1);


Coef=zeros(NumA+NumB,Deg+1);
M=7;
Fib=[1,0;1,1];
for i=3:(M+Deg)
    [rows, cols]=size(Fib);
    z1=zeros(rows,1);
    Fib=[Fib,z1];
    last=Fib(end,:);
    nr=zeros(1, cols+1);
    nr(1)=1;
    for i=1:cols
        nr(1+i)=last(i)+last(i+1);
    end
    Fib=[Fib;nr];
end
Fib 
%listSStar=zeros(1,NumPeriods);
%tic;
InitFFT_back(modelType);
%t=toc;
%sprintf('time spent for InitFFT : %f',t)
%exit(0)
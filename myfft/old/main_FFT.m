

%function  main_FFT(model, scenario)
function  main_FFT()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

model = 'Gauss';
scenario = 1;
input_path=strcat('');
input_fileName='input.xlsx'; 

[a, b] =xlsread(strcat(input_path,input_fileName),model);
[l,w]=size(a);
if scenario > l
    error('invalid scenario %d', 2);
end
inputVar = containers.Map(b,a(scenario,:));
inputVar = containers.Map(b,a(scenario,:));
NInter = inputVar('NInter');
TT = inputVar('TT')/NInter;
Delta = TT/NInter;
Rate = inputVar('r');
Sigma = inputVar('sigma');
Discount=exp(-Rate*Delta);
[values, index] =xlsread(strcat(input_path,input_fileName),'Strikes');
S = values(1);
K = values(2:end);

[values, index] =xlsread(strcat(input_path,input_fileName),'FFT'); 
FFTVar = containers.Map(index,values(scenario,:));
fftpower = FFTVar('fftpower');
Left = FFTVar('Left');
Right = FFTVar('Right'); 
M_A = FFTVar('M_A');
M_D = FFTVar('M_D');
Oalpha = FFTVar('alpha');
Obeta = FFTVar('ind');

N=2^fftpower;
F=DensityCoeficient(model,inputVar, Left,Right, N);
[D_fft, AD_fft] = FFT(F,N, M_A, M_D, Left, Right);
Exp_fft = FFTExp(1,F, N, Left, Right);
x = MyChar(model,inputVar,-complex(0,1)*1) - Exp_fft(N/2);

Deg = 2;

Lambda=(Right-Left)/(N/2); %new
D=(0:1:N/2)*Lambda+Left;
Euro=zeros(1,N/2+1);
CP=N/4+1; %center point index
%D(k) = (k-1)*Lambda + Left for k = 1, 2, ..., N/2
%D(CP)= (Left+Right)/2.   if Left + Right = 0, then D(CP) = 0
NumA = 3; %how many piecewise functions are used for a interval ending at D(CP)
SA = 2; % the grid points we choose for interpolation.  SA = 2, we select every other points
st=CP-Deg*NumA*SA; %starting point position
GP=st+(0:1:Deg*NumA)*SA; % poistions of the points that are used for interpolation
B=st+(0:1:NumA)*SA*Deg; % the boundary points for those intervals over which a polynomial is setup;

NumB=3; %similar to NumA, but for the interval [D(CP), Right]
SB=floor((N/4)/(Deg*NumB));
gbp=CP+(1:1:Deg*NumB)*SB;
GP=[GP,gbp];
bb=CP+(1:1:NumB)*Deg*SB;
B=[B,bb];
binCoef = Fib(M_D + Deg);

test_flag_Euro = 1;
%% Test Euro
if test_flag_Euro == 1
    low=min(max(D(B(1))-D(GP), D(1)), D(N/2));
    lowpos=min(max(B(1)-GP+CP,1),N/2);
    Strike = S;
    Euro2=Discount*(AD_fft(1,lowpos)*Strike*Oalpha-Obeta*exp(D(GP)).*Exp_fft(lowpos)); %new

    low=min(max(D(CP)-D(GP), D(1)), D(N/2));
    lowpos=min(max(CP-GP+CP,1),N/2);
    Euro=Discount*(AD_fft(1,lowpos)*Strike*Oalpha-Obeta*exp(D(GP)).*Exp_fft(lowpos)); %new

    [Euro2', Euro']


    x=D(GP);
    val=BS(x, K, Rate, Delta, Sigma);
    max(abs(val-Euro))
end

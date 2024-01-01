%test PIntA and PIntB

function void =testPartIntV2()
global Sigma Rate Delta D N Left Right Lambda PIntA PIntB Deg SA SB%for testing purpose
mu=(Rate-0.5*Sigma^2)*Delta;
std=Delta^0.5*Sigma;
%delta=delta*Lambda;


%lambdaA=Deg*SA*Lambda;
%lambdaB=Deg*SB*Lambda;
poslow=(1:1:N/2);
posupA=poslow+Deg*SA;
posupA=min(posupA,N/2);
posupB=poslow+Deg*SB;
posupB=min(posupB,N/2);

%up=D(2:N/2);
%d=D(1:N/2-1);%up;%-Lambda/2;
%low=d+delta;

up=D(posupB);
low=D(1:N/2);

%delta=0;%Lambda*0.5;

%delta=Lambda*SA*Deg/2
disp('delta is')
delta=Lambda*SB*Deg
up=(up-mu)/std;
%if delta==0, then testing on grid points
low=(low+delta-mu)/std;
me=mu-D(1:N/2); %mean



val1=std*GetNormPolyInt(up,low,1)+me.*GetNormPolyInt(up,low,0);
%if delta==0
%    [PIntA(1,1:N/4)',val1(1:N/4)']
%    max(abs(PIntA(1,:)-val1))
%else
    testv=zeros(1,N/2);
    for i=1:N/2
        testv(i)=GetPartInt(delta, i, 1, 2);
    end
    %[testv(1,1:N/4)', val1(1:N/4)']
    max(abs(testv-val1))
%end
val2=std^2*GetNormPolyInt(up,low,2)+...
    2*std*me.*GetNormPolyInt(up,low,1)+...
    me.^2.*GetNormPolyInt(up,low,0);
%if delta==0
%    max(abs(PIntA(2,:)-val2))
%else
    testv=zeros(1,N/2);
    for i=1:N/2
        testv(i)=GetPartInt(delta, i, 2, 2);
    end
    max(abs(testv-val2))
%end

val3=std^3*GetNormPolyInt(up ,low , 3)+3*std^2*me.*GetNormPolyInt(up,low,2)+...
    3*std*me.^2.*GetNormPolyInt(up,low,1)+me.^3.*GetNormPolyInt(up,low,0);
%if delta==0
%    max(abs(PIntA(3,:)-val3))
%else
    testv=zeros(1,N/2);
    for i=1:N/2
        testv(i)=GetPartInt(delta, i, 3, 2);
    end
    max(abs(testv-val3))
%end


val4=std^4*GetNormPolyInt(up ,low , 4)+4*std^3*me.*GetNormPolyInt(up,low,3)+...
                6*std^2*me.^2.*GetNormPolyInt(up,low,2)+4*std*me.^3.*GetNormPolyInt(up,low,1)+...
                me.^4.*GetNormPolyInt(up,low,0);            
%if delta==0
%    max(abs(PIntA(4,:)-val4))
%else
    testv=zeros(1,N/2);
    for i=1:N/2
        testv(i)=GetPartInt(delta, i, 4, 2);
    end
    max(abs(testv-val4))
%end

%[SA,SB]
%test results for delta =0, numperiods=10 left=-1, right=1  M=7, Deg=4
%fftpower=7  [7.10e-017, 1.94e-016, 4.16e-016 5.83e-016]
%fftpower=8  [1.04e-016, 3.05e-016, 2.98e-016, 4.44e-016]

%test on nongrid periods on SA

%delta =0.03125=Deg*SA*Lambda fftpower=9
%[1.05e-010, 2.96e-012, 8.41e-014, 2.51e-015]

%for fftpower=8, delta=0.0625
%[5.36e-008, 3.01e-009, 1.71e-010, 9.80e-012]

%fftpower=10, delta=0.015625, 
%[2.06e-013, 2.91e-015, 6.24e-016, 8.75e-016]

%fftpower=11, 0.0078125,
%[3.31e-015, 4.12e-016, 5.74e-016, 1.13e-015

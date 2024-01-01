function void =testPartInt();
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
low=D;

x=zeros(1,N/2);
y=zeros(1,N/2);

up=(up-mu)/std;
low=(low-mu)/std;
me=mu-D;

for deg=1:4
    for i=1:N/2
        %if deg==0
        %    x(i)=GetNormPolyInt((up(i)-mu)/std,(low(i)-mu)/std,0);
        %    x(i)=GetDerivative(up(i), 0)-GetDerivative(low(i), 0);
        if deg==1
            x(i)=std*GetNormPolyInt((up(i)-mu)/std,(low(i)-mu)/std,1)+...
                (mu-low(i))*GetNormPolyInt((up(i)-mu)/std,(low(i)-mu)/std,0);
        elseif deg==2
            x(i)=std^2*GetNormPolyInt((up(i)-mu)/std,(low(i)-mu)/std,2)+...
                2*std*(mu-low(i))*GetNormPolyInt((up(i)-mu)/std,(low(i)-mu)/std,1)+...
                (mu-low(i))^2*GetNormPolyInt((up(i)-mu)/std,(low(i)-mu)/std,0);
        elseif deg==3
            val3=std^3*GetNormPolyInt(up ,low , 3)+3*std^2*me.*GetNormPolyInt(up,low,2)+...
                3*std*me.^2.*GetNormPolyInt(up,low,1)+me.^3.*GetNormPolyInt(up,low,0);
            max(abs(pint(4,1:N/2-1)-val3))
        elseif deg==4
            val4=std^4*GetNormPolyInt(up ,low , 4)+4*std^3*me.*GetNormPolyInt(up,low,3)+...
                6*std^2*me.^2.*GetNormPolyInt(up,low,2)+4*std*me.^3.*GetNormPolyInt(up,low,1)+...
                me.^4.*GetNormPolyInt(up,low,0);
            max(abs(pint(5,1:N/2-1)-val4))
        end
        y(i)=PIntB(deg, i);
    end
    %disp('really')
    %[x',y']
    max(abs(x-y))
end
%[SA,SB]
%test results 
%fftpwer=5,  [0.00203150185136  9.828583204476782e-004]
%fftpower=6  [6.621723646627054e-007, 1.617437174533318e-007]
%fftpower=7  [2.636779683484747e-016, 3.053113317719181e-016]

%fftpower=7, [2.636779683484747e-016, 3.053113317719181e-016]
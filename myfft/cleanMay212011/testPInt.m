function void =testPInt(model)
global Sigma Rate Delta D N Left Right  %for testing purpose
mu=(Rate-0.5*Sigma^2)*Delta;
std=Delta^0.5*Sigma;

up=(D(2:N/2)-mu)/std;
low=(D(1:N/2-1)-mu)/std;
me=mu-D(1:N/2-1);


%fftad=GetAntiDerivative(model);
pint=GetPartInt(model);

val0=GetNormPolyInt(up ,low , 0);
max(abs(pint(1,1:N/2-1)-val0))

val1=std*GetNormPolyInt(up,low,1)+me.*GetNormPolyInt(up ,low , 0);
max(abs(pint(2,1:N/2-1)-val1))

val2=std^2*GetNormPolyInt(up ,low , 2)+2*std*me.*GetNormPolyInt(up,low,1)+...
    me.^2.*GetNormPolyInt(up,low,0);
max(abs(pint(3,1:N/2-1)-val2))

val3=std^3*GetNormPolyInt(up ,low , 3)+3*std^2*me.*GetNormPolyInt(up,low,2)+...
    3*std*me.^2.*GetNormPolyInt(up,low,1)+me.^3.*GetNormPolyInt(up,low,0);
max(abs(pint(4,1:N/2-1)-val3))

val4=std^4*GetNormPolyInt(up ,low , 4)+4*std^3*me.*GetNormPolyInt(up,low,3)+...
    6*std^2*me.^2.*GetNormPolyInt(up,low,2)+4*std*me.^3.*GetNormPolyInt(up,low,1)+...
    me.^4.*GetNormPolyInt(up,low,0);
max(abs(pint(5,1:N/2-1)-val4))
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
%exit('testDensity')
function void =testDensity(numPeriods)
global Sigma Rate Delta D N  M  Left Right  %for testing purpose
%main(1, 90, 8, 4,7, -1, 1, numPeriods)
Constant(1,numPeriods);
Init(1, 90, 8, 4,7, -1, 1, numPeriods) 
K=0:1:N-1;
F=DensityCoeficient(model, Left,Right, N); 
mu=(Rate-0.5*Sigma^2)*Delta;
std=Delta^0.5*Sigma;
val0=exp(-0.5*((D-mu)/std).^2)/((2*pi)^0.5*std);
val1=zeros(1,N/2+1);
for i=1:N/2+1;
    val1(i)=sum(F.*cos(K*pi*(D(i)-Left)/(Right-Left)));
end
max(abs(val0-val1))
scaling=std;
F=F/ModelScaling; %this will be FFT for Y=ModelScaling *X
[max(F), min(F)]
exit('temp in Init')
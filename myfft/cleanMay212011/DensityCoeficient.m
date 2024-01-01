function coefdensity=DensityCoeficient(modelType, leftend,rightend, N)
x=0:1:(N-1);
coefdensity=MyChar(modelType, x*pi/(rightend-leftend)).*exp(-x*leftend*pi/(rightend-leftend)*complex(0,1));
coefdensity=2*real (coefdensity)/(rightend-leftend);
coefdensity(1)=coefdensity(1)/2;
function coefdensity=DensityCoeficient(model, inputVar, leftend,rightend, N)
%Formula (7) on page 2.  where
% l = rightend-leftend
% L = leftend
% coefdensity = {\bar{A}_j}_{j=0,...,N-1}
x=0:1:(N-1);
coefdensity=MyChar(model, inputVar,  x*pi/(rightend-leftend)).*exp(-x*leftend*pi/(rightend-leftend)*complex(0,1));
coefdensity=2*real (coefdensity)/(rightend-leftend);
coefdensity(1)=coefdensity(1)/2;
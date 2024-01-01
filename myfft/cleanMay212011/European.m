function void =European(modelType)
global D Euro Discount Oalpha Obeta Strike N FFTExp FFTAD Euro2 GP B CP


lowpos=min(max(B(1)-GP+CP,1),N/2);
Euro2=Discount*(FFTAD(1,lowpos)*Strike*Oalpha-Obeta*exp(D(GP)).*FFTExp(lowpos)); %new



lowpos=min(max(CP-GP+CP,1),N/2);
Euro=Discount*(FFTAD(1,lowpos)*Strike*Oalpha-Obeta*exp(D(GP)).*FFTExp(lowpos)); %new



x=D(GP);
val=BS(x);
max(abs(val-Euro))
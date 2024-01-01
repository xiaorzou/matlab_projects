function val=Continuation()
global Oalpha Obeta Strike Discount NumA NumB Coef 
global NumPeriods EpsilonSearching Lambda B N 
global CP D Deg FFTAD SA Euro Euro2 FFTExp D3
global Index GP Left Right %donot need this one late.

if Index==NumPeriods-1
    val=Euro;%max(Euro(GP),0);
    return;
end

if Index==-1
    disp('in conti')
val(1)
end

[row,cols]=size(GP);
val1=0;
for pos=1:cols
    val1(pos)=sum(sum(Coef(:,(1:Deg)).*D3(:,:,pos)));
    %[Coef(1:NumA,(1:Deg)), D3(1:NumA,:,pos)]
    %exit(0)
    lowpos=min(max(B(1:NumA+NumB)-GP(pos)+CP,1), N/2);
    uppos=min(max(B(2:NumA+NumB+1)-GP(pos)+CP,1),N/2);
    %sum(Coef(:,Deg+1)'.*(FFTAD(1,uppos)-FFTAD(1,lowpos)))
    %exit('dd')
    val1(pos)=val1(pos)+sum(Coef(:,Deg+1)'.*(FFTAD(1,uppos)-FFTAD(1,lowpos)));
end
val1=max(val1,0);
%[val1(1),Euro2(1)]
val=Euro2+val1*Discount;
%val(1)
%exit(0)
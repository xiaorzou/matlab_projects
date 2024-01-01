function void =Update()  %timeIndex=NumPeriods-1:-1:0
global NumPeriods Coef NumA NumB N CP D Deg Lambda
global Index GP B CP SA

coef=zeros(NumA+NumB,Deg+1);
payoff=zeros(1,Deg+1);
flag=0;
cont=Continuation();
cc=cont(1);
pp=GetPayOff(D(B(1)));

if cc>pp
    error('cont(1)=%f should smaller than payoff=%f, need smaller SP. Index=%f', cc, pp, Index);
end

for i=1:(NumA+NumB);
    pos=GP((1+(i-1)*Deg):(1+i*Deg));
    val=max(cont((1+(i-1)*Deg):(1+i*Deg)), GetPayOff(D(pos)));
    coef(i,:)=GetCoef(D(pos),val, Deg);
end
if Index==-1 %check
    for i=1:(NumA+NumB)
        [max(abs(polyval(coef(i,:), D(GP((1+(i-1)*Deg):(1+i*Deg)))-D(GP(1+(i-1)*Deg)))-...
            BS(D(GP((1+(i-1)*Deg):(1+i*Deg)))))),max(abs(polyval(coef(i,:), D(GP((1+(i-1)*Deg):(1+i*Deg)))-D(GP(1+(i-1)*Deg)))-...
            GetPayOff(D(GP((1+(i-1)*Deg):(1+i*Deg))))))]
        if i==8
            cont((1+(i-1)*Deg):(1+i*Deg))
            polyval(coef(i,:), D(GP((1+(i-1)*Deg):(1+i*Deg)))-D(GP(1+(i-1)*Deg)))
            GetPayOff(D(GP((1+(i-1)*Deg):(1+i*Deg))))
            %exit('kkkkkkkk')
        end 
    end
    %exit('tempmmm')
end

%update information here, not before
Coef=coef;
%exit(0)
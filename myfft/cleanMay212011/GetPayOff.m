%function val=GetPayOff(pos)
function val=GetPayOff(logPrice)
global Oalpha Obeta Strike 
s=exp(logPrice);
val=0;
if s<=Strike
    val=max(Oalpha*Strike-Obeta*s,0);
end

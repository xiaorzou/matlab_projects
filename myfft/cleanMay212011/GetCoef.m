function coefb=GetCoef(x,cont, deg)

%test part
%coef=[1,-1,1,-1,1];
%x=[1,2,3,4,5]
%y=polyval(coe,x,4)
%deg=4;
%end of test part


coef=polyfit(x,cont,deg);

%code
g=1;
coefb=zeros(1,deg+1);
for i=(deg+1):-1:1
    if i~=deg+1
        g=g*(deg+1-i);
    end
    %[i-1, g]
    %coef
    coefb(i)=polyval(coef,x(1), i-1)/g;
    coef=coef.*(i-1:-1:0);
    coef=coef(1:i-1);
end
%coefb
%polyval(coefb,x-x(1),deg)
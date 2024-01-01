%the polynomial is coef(1)(x-left)^deg+...
function x=Newton(coef, left, right)
global Oalpha Obeta  Strike Deg EpsilonSearching

k=Deg:-1:0;
coe=coef.*k;
coe=coe(1:Deg);
expb=exp(left);
if (polyval(coef,0)-GetPayOff(left)>0) | (polyval(coef,right-left)-GetPayOff(right)<0)
    coef
    bpos
    error('wrong in Newton')
end
x=left+(right-left)/2;
y=polyval(coef,x-left)-GetPayOff(x);
%y=polyval(Coef,x);
counter=1;
while abs(y)>EpsilonSearching
    x=x-y/(polyval(coe,x-left)+Obeta*exp(x));
    y=polyval(coef,x-left)-Oalpha*Strike+Obeta*exp(x);
    counter=counter+1;
    %y
    if counter>10
        [y]
        error('warning in newton method, searching more than 10 times');
    end
end
if x<left | x>right
    disp('wrong in Netwon, x, guessleft, guessright are');
    [x, left, right]
    [polyval(coef, x-left), Oalpha*Strike-Obeta*exp(x)]
    [polyval(coef,0), Oalpha*Strike-Obeta*exp(left), polyval(coef, right-left), Oalpha*Strike-Obeta*exp(right)] 
    
    
    exit(1)
end
%counter
%x
%[polyval(coef, x-left), BS(x)]
%exit('in newton')
%polyval(coef,x-left)-GetPayOff(x)
%exit('in newton')
%error('temp in newton')

%x=-0.02501697449296 for 2 period


    



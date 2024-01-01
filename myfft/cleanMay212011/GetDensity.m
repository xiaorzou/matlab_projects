%q(x)^(order)
function val=GetDensity(x, order) %order-1 is the degree of derivative, order=1, ..., M-1
global D FFTF Lambda M Left Right N

if order>=M
    error('in GetDensity, the order %d is bigger than M-1:%d ', order, M-1);
end
if (x<Left | x>Right) %& order==0 , for order >0, if x is out of range, what we should return?
    val=0;
    return;
end
q=floor((x-Left)/Lambda);
r=x-Left - Lambda *q;
if q>=N/2
    q=N/2-1;
    r=0;
end
if (r/Lambda)<10^(-14)
%    disp('enter Gend1');
    val=FFTF(order,q+1);
    %[order,val]
    return;
elseif (r/Lambda)>(1-10^(-14))
    val=FFTF(order, q+2);
    return;
else
    pos=0;
    if (r/Lambda)<0.5
        pos=q+1;
    else
        pos=q+2;
    end
    delta=x-D(pos);
    c=1;
    val=FFTF(order,pos);
    for i=1:M-order
        c=c*delta/i;
        val=val+c*FFTF(i+order,pos);
    end
end


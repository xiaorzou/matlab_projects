%  order=0, 1, 2, 3, ... M
%  GetDerivative(x, 0)=CDF(x)
%  GetDerivative(x, 1)=q(x)
%  GetDerivative(x, order)=GetDerivative(x, order-1)'
function val=GetDerivative(x, order)
global D FFTF Lambda M Left Right N FloatingError CDF FFTAD

if order>M
    error('in GetDensity, the order %d is bigger than M:%d ', order, M);
end
if (x<Left | x>Right)
    if order==1
        val=0; %density
    elseif order==0 %CDF
        if x>Right 
            val=1;
        elseif x<Left
            val=0;
        end
    else %for other order, we are not sure how to set the value.
        error('wrong range for x in GetDerivative');
    end
    return;
end
q=floor((x-Left)/Lambda);
r=x-Left - Lambda *q;
if q>=N/2
    q=N/2-1;
    r=0;
end
if (r/Lambda)<FloatingError
%    disp('enter Gend1');
    if order==0
        %val=CDF(q+1);
        val=FFTAD(1,q+1);
    else
        val=FFTF(order,q+1);
    end
    return;
elseif (r/Lambda)>(1-FloatingError)
    if order==0
        %val=CDF(q+2);
        val=FFTAD(1,q+2);
    else
        val=FFTF(order, q+2);
    end
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
    if order==0 %for cdf
        %val=CDF(pos);  
        val=FFTAD(1,pos);
        %+delta*FFTF(1,lpos)+delta^2*FFTF(2,lpos)/2+delta^3*FFTF(3,lpos)/6+delta^4*FFTF(4, lpos)/24 +...
        %delta^5*FFTF(5, lpos)/120+delta^6*FFTF(6, lpos)/720+delta^7*FFTF(7, lpos)/5040;
    else
        val=FFTF(order,pos);
    end
    for i=1:M-order
    %for i=1:6
        c=c*delta/i;
        val=val+c*FFTF(i+order,pos);
    end
end

%get CDF
%if pos>0, then CDF at D(pos)+offset
%if pos<=0, then CDF at x

function val=GetCDF(x, offset,pos)
global D FFTF FFTAD Lambda Left Right N FloatingError

if pos<1
    if x<Left | x>Right
        error('Wrong range for x %d', x)
    end
    q=floor((x-Left)/Lambda);
    r=x-Left - Lambda *q;
    if q==N/2
        q=q-1;
    end
    lpos=q+1;
    if r> (Lambda/2)
        lpos=lpos+1;
    end
    delta=x-D(lpos);
    if abs(delta)>Lambda
        error('in GetAntiDerivative')
    end
else
    delta=offset;
    lpos=pos;
end
if abs(delta)<FloatingError
    val=FFTAD(1,lpos);
else
    val=FFTAD(1,lpos)+delta*FFTF(1,lpos)+delta^2*FFTF(2,lpos)/2+delta^3*FFTF(3,lpos)/6+delta^4*FFTF(4, lpos)/24 +...
        delta^5*FFTF(5, lpos)/120+delta^6*FFTF(6, lpos)/720+delta^7*FFTF(7, lpos)/5040;
end
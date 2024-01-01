%q_{(order)}(at): if order=0, the value of q(at), otherwise the value of
%antiderivative at 'at', expend around D(lpos)
%order =0, 1, 2, ..., 5 
%q(at)=GetAntiDerivative(at, lpos, 0)
%\int^at_{Left} q(x)dx=GetAntiDerivative(at, lpos, 1)

function val=GetAntiDerivative(x, order)
global D FFTF FFTAD Lambda

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
if order==0
    if abs(delta)<FloatingError
        val=FFTF(1, lpos);
    else
        val=FFTF(1, lpos)+delta*FFTF(2,lpos)+0.5*delta^2*FFTF(3,lpos)+delta^3*FFTF(4,lpos)/6+...
            delta^4*FFTF(5, lpos)/24 +delta^5*FFTF(6, lpos)/120 +delta^6*FFTF(7, lpos)/720;
    end
elseif order ==1 % \int^at_{-\infty} q(x)dx
    if abs(delta)<FloatingError
        val=FFTAD(1,lpos);
    else
        val=FFTAD(1,lpos)+delta*FFTF(1,lpos)+delta^2*FFTF(2,lpos)/2+delta^3*FFTF(3,lpos)/6+delta^4*FFTF(4, lpos)/24 +...
        delta^5*FFTF(5, lpos)/120+delta^6*FFTF(6, lpos)/720+delta^7*FFTF(7, lpos)/5040;
    end
elseif order==2
    if abs(delta)<FloatingError
        val=FFTAD(2,lpos);
    else
        val=FFTAD(2,lpos)+delta*FFTAF(1,lpos)+delta^2*FFTF(1,lpos)/2+delta^3*FFTF(2,lpos)/6 + ...
        delta^4*FFTF(3, lpos)/24 +delta^5*FFTF(4, lpos)/120 +...
        delta^6*FFTF(5, lpos)/720+delta^7*FFTF(6, lpos)/5040+delta^8*FFTF(7, lpos)/40320;
    end
elseif order==3
    if abs(delta)<FloatingError
        val=FFTAD(3,lpos);
    else
        val=FFTAD(3,lpos)+delta*FFTAD(2,lpos)+delta^2*FFTAD(1,lpos)/2+delta^3*FFTF(1,lpos)/6 + ...
        delta^4*FFTF(2, lpos)/24 +delta^5*FFTF(3, lpos)/120 +...
        delta^6*FFTF(4, lpos)/720+delta^7*FFTF(5, lpos)/5040+delta^8*FFTF(6, lpos)/40320+...
        delta^9*FFTF(7, lpos)/362880;
    end
elseif order==4
    if abs(delta)<FloatingError
        val=FFTAD(4,lpos);
    else
        val=FFTAD(4,lpos)+delta*FFTAD(3,lpos)+delta^2*FFTAD(2,lpos)/2+delta^3*FFTAD(1,lpos)/6 + ...
        delta^4*FFTF(1, lpos)/24 +delta^5*FFTF(2, lpos)/120 +...
        delta^6*FFTF(3, lpos)/720+delta^7*FFTF(4, lpos)/5040+delta^8*FFTF(5, lpos)/40320+...
        delta^9*FFTF(6, lpos)/362880+delta^10*FFTF(7, lpos)/3628800;
    end
elseif order==5
    if abs(delta)<FloatingError
        val=FFTAD(5,lpos);
    else
        val=FFTAD(5,lpos)+delta*FFTAD(4,lpos)+delta^2*FFTAD(3,lpos)/2+delta^3*FFTAD(2,lpos)/6 + ...
        delta^4*FFTAD(1, lpos)/24 +delta^5*FFTF(1, lpos)/120 +...
        delta^6*FFTF(2, lpos)/720+delta^7*FFTF(3, lpos)/5040+delta^8*FFTF(4, lpos)/40320+...
        delta^9*FFTF(5, lpos)/362880+ delta^10*FFTF(6, lpos)/3628800 +delta^11*FFTF(7, lpos)/39916800;
    end
else
    error('not implemented yet, anti derivative order up to 5 in GetAntiDerivative');
end
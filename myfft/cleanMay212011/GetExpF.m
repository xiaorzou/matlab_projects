%int^x_Left e^(hu) q(u)du   for x\in [Left Right]
%only implement h==1!
function val = GetExpF(x, h)
global N Left Right Lambda  FFTF FFTExp D M  FloatingError
if x<Left | x>Right
    error('Wrong range for x %d', x)
end
if h~=1
    error('has not implemented for degree bigger than 1. the degree is %f', h);
end
q=floor((x-Left)/Lambda);
r=x-Left - Lambda *q;
if q==N/2
    q=q-1;
end
if abs(r/Lambda)<FloatingError
    val=FFTExp(q+1);
    return;
end

if abs(r/Lambda)>(1-FloatingError)
    val=FFTExp(q+2);
    return;
end

pos=q+1;
if r> (Lambda/2)
    pos=pos+1;
end
del=x-D(pos);
val=FFTExp(pos);
c=1;
for i=1:M-1
    c=c*del/i;
    val=val+c*getDerAtGrid(pos,i);
end

function val=getDerAtGrid(pos,n)
global D FFTF
val=0;
for j=1:n
    val=val+FFTF(j, pos)* nchoosek(n-1,j-1);
end
val=val*exp(D(pos));
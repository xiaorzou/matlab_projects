% \int^(D(base+del)-D(pos))_(y+D(pase)-D(pos)) (x-(D(base)-D(pos)))^deg q(x) 
%where del=SA*Deg if option ==1 or SB*Deg if option==2
%deg=0,1,2,3,...,Deg, pos=1,2,...,N/2, base points are on GPA or GPB
function val=GetPartInt(y, base, pos, deg, option) 
global N Deg D Left Right Lambda FFTF  PIntA PIntB FloatingError SA SB M Fib CP


%disp('enter !!!')
if deg<0 | deg>Deg
    error('wrong deg value in GetPartInt, %f', deg);
end
del=0;
if option == 1 %& 
    del=SA*Deg;
    if y<0   | y >Lambda*SA*Deg
        error('wrong y value in GetPartInt, %f', y);
    end
elseif option == 2
    del=SB*Deg;
    if y<0   | y >Lambda*SB*Deg
        error('wrong y value in GetPartInt, %f', y);
    end
else
   error('wrong option value in GetPartInt, %f', option);
end

posNew=base-pos+CP;
posup =max(min(posNew+del, N/2),1);
poslow=max(min(posNew,N/2),1);

%disp('31')
if deg==0
    val=FFTAD(1,posup)-FFTAD(1,posup);
else
    val=GetPIntNew(base, pos, deg, option);
end
if abs(y)<FloatingError | posNew<1  %if posNew<1, then all the derivatives of q(x) are 0 at D(1)!
    return;
end
%deg=1
%M=7;

%if deg==1
%    val=val-0*y-FFTF(1,pos)*y^2/2-2*FFTF(2,pos)*y^3/6-3*FFTF(3,pos)*y^4/24-...
%        4*FFTF(4,pos)*y^5/120-5*FFTF(5,pos)*y^6/720-6*FFTF(6,pos)*y^7/40320-...
%        7*FFTF(7,pos)*y^8/362880;
%else 

%myfact=1;
myfact=GetFactor(deg);
val=val-myfact*FFTF(1, posNew)*y^(deg+1)/(myfact*(deg+1));
g=myfact*(deg+1);
%c=myfact;
%Fib =[1     0     0     0     0     0     0     0     0     0     0
%     1     1     0     0     0     0     0     0     0     0     0
%     1     2     1     0     0     0     0     0     0     0     0
%     1     3     3     1     0     0     0     0     0     0     0
%     1     4     6     4     1     0     0     0     0     0     0
%     1     5    10    10     5     1     0     0     0     0     0
%     1     6    15    20    15     6     1     0     0     0     0
%     1     7    21    35    35    21     7     1     0     0     0
%     1     8    28    56    70    56    28     8     1     0     0
%     1     9    36    84   126   126    84    36     9     1     0
%     1    10    45   120   210   252   210   120    45    10     1];
for n=(deg+2):(M+deg)
    g=g*n;
    val=val-Fib(n,deg+1)*myfact*FFTF(n-deg, posNew)*y^(n)/g;
    %c=[c,Fib(n,deg+1)*myfact];
    %[n,g]
end
%c
%end




function val=GetFactor(n)
if n==1
    val=1;
elseif n==2
    val=2;
elseif n==3
    val=6;
elseif n==4
    val=24;
elseif n==5
    val=120;
elseif n==6
    val=720;
elseif n==7
    val=5040;
elseif n==8
    val=40320;
elseif n==9
    val=362880;
elseif n==10
    val=3628800;
else
    val=n*GetFactor(n-1);
end
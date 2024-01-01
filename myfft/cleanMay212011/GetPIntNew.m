%pos can be a vector
%val=\int^(D(base+Deg*SA)-D(pos))_(D(base)-D(pos))(x+D(pos)-D(base))^deg q(x)dx
%treat it as 
%\int^(upPos)_(lowPos)(x+D(pos)-D(base))^degq(x)dx     
%where upPos=min(max (base+Deg*SA-pos+CP, 1), N/2);
%      lowPos=min(max (base-pos+CP, 1), N/2);
function val=GetPIntNew(base, pos, deg, option)

global CP D Deg SA SB Lambda FFTAD PIntA PIntB N

posNew=base-pos+CP;
if (posNew>0) & (posNew<N/2+1)
    if option==1
        val=PIntA(deg, posNew);
    elseif option==2
        val=PIntB(deg, posNew);
    else
         error('wrong1 option in GetPIntNew');
    end
else
    del=Deg*SA;
    if option==2
        del=Deg*SB;
    elseif option ~= 1
        error('wrong2 option in GetPIntNew');
    end
    upPos=min(max(posNew+del,1),N/2);
    lowPos=min(max(posNew,1),N/2);
    val=FFTAD(1,upPos)*(D(upPos)+D(pos)-D(base))^deg - FFTAD(1,lowPos)*(D(lowPos)+D(pos)-D(base))^deg;
    factor=1;
    for k=1:deg
        factor=factor*(deg-k+1);
        val=val+(-1)^k*factor*( FFTAD(k+1,upPos).*( D(upPos)+D(pos)-D(base)).^(deg-k)-...  %lambdaA^(deg-k)-...
            FFTAD(k+1,lowPos)*(D(lowPos)+D(pos)-D(base))^(deg-k));
    end
end

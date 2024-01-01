%
%function val=GetPIntA(Deg-h+1, BA(i)-pos+CP)

%\int^(D(pos+Deg*SA))_(D(pos))(x-D(pos))^deg q(x)dx
%need carefully handle  pos when pos<0 and pos+Deg*SA>0
%\int^(D(base+Deg*SA)-D(pos))_(D(base)-D(pos))(x-(D(base)-D(pos)))^deg q(x)dx
function val=GetPInt(deg, base, pos, option)
global CP D Deg SA Lambda FFTAD PIntA 
if option==1
    del=Deg*SA;
elseif option==2
    del=Deg*SB;
else
    error('wrong option in GetPInt');
end

if D(base)-D(pos)>=D(1)
    if option==1
        val=PIntA(deg, base-pos+CP);
    else %option==2
        val=PIntB(deg, base-pos+CP);
    end
else %pos<0, 
    upValue=min(max(D(base+del)-D(pos),D(1)),D(N/2));
    lowValue=min(max(D(base)-D(pos),D(1)),D(N/2));
    upPos=min(max(base+del-pos+CP,1),N/2);
    lowPos=min(max(base-pos+CP,1),N/2);
        %poslow=1;
        %posupA=pos+del;
    val=FFTAD(1,upPos)*(D(upPos)+D(pos)-D(base))^deg - FFTAD(1,lowPos)*(D(lowPos)+D(pos)-D(base))^deg;
    factor=1;
    for k=1:deg
        factor=factor*(deg-k+1);
        val=val+(-1)^k*factor*( FFTAD(k+1,upPos).*( D(upPos)+D(pos)-D(base)).^(deg-k)-...  %lambdaA^(deg-k)-...
            FFTAD(k+1,lowPos)*(D(lowPos)+D(pos)-D(base))^(deg-k));
    end
end

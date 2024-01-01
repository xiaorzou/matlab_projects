function val=GetNormPolyInt(up,low,deg) %\int^up_low x^deg e^{-x^2/2}/(2pi)^0.5 dx

if deg==0
    val=normcdf(up)-normcdf(low);
elseif deg==1
    val=-(exp(-0.5*up.^2)-exp(-0.5*low.^2))/(2*pi)^0.5;
else
    val=-(up.^(deg-1).*exp(-0.5*up.^2)-low.^(deg-1).*exp(-0.5*low.^2))/(2*pi)^0.5+...
        (deg-1)*GetNormPolyInt(up, low, deg-2);
end
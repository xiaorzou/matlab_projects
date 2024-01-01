%s=exp(x) is the stock price!

%function vector=BS(x)

function val=BS(x)
global Delta Rate Sigma Strike
s=exp(x);
vector=Call(s);
val=vector-s+Strike*exp(-Rate*Delta);
%vector(1)=vector(1)-s+Strike*exp(-Rate*Delta);
%vector(2)=vector(2)-1;


function vector=Call(s)
global Delta Rate Sigma Strike
d1=(log(s/Strike)+(Rate+Sigma^2/2)*Delta)/(Sigma*Delta^0.5);
d2=d1-Sigma*Delta^0.5;
val=s.*normcdf(d1)-Strike*exp(-Rate*Delta).*normcdf(d2);
%delta=normcdf(d1);
%gamma=exp(-0.5*d1.^2)/((2*pi*Delta)^0.5*Sigma*s);
%vector=[val, delta, gamma];
vector=val;
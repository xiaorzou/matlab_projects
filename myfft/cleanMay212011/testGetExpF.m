%int^x_Left e^(h*u) q(u)du   for x\in [Left Right]

%test 
function void = testGetExpF(h)
global Sigma Rate Delta D N  M  Left Right  %for testing purpose
mu=(Rate-0.5*Sigma^2)*Delta
std=Delta^0.5*Sigma
xvalue=(D(1:(N/2-1))+D(2:(N/2)))/2;
%xvalue=D(1:(N/2-1));
yvalue=zeros(1,N/2-1);
for i=1:(N/2-1)
    yvalue(i)=GetExpF(xvalue(i),h);
end
v=exp(h*mu+0.5*std^2*h^2)*(normcdf((xvalue-mu)/std-h*std)-normcdf((Left-mu)/std-h*std));
max(abs(v-yvalue))
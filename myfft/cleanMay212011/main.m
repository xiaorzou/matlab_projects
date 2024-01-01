%function void = main(model, strike, fftpower,  skip1,skip2, deg,m, left, right)

function val = main(model, strike, fftpower, deg,m, left, right, numPeriods,alpha, beta, numA, numB)

global KStar NumPeriods ScaleUnit Strike Lambda Coef Deg  NumA NumB Index B
global bug M D% for debug, remove late 

%T=0
%for i=1:10000
%tic;
%Init(model, strike, fftpower, skip1, skip2,  deg,m, left, right);
Constant(model,numPeriods, alpha, beta);
Init(model, strike, fftpower, deg,m, left, right, numA, numB);
European();
testDerivative();
%testPartInt();
%testPInt(model);
testDensity();
h = 1;
testGetExpF(h)

for Index=NumPeriods-1:-1:0
    bug=NumPeriods-2;
    Update();
end
%exit('temp exit in main');
%x=[36,38,40,42,44]';
x=100;
[rows, cols]=size(x);
y=zeros(rows,1);
c=zeros(rows,1);
logStock=log(x/ScaleUnit);
%exp(logStock)*ScaleUnit
counter=1;
%Coef
for k=1:rows
    if logStock(k)<=D(B(1))
        c(k)=GetPayOff(logStock(k));
    else
        flag=0;
        for m=1:NumA+NumB
            if D(B(m))<=logStock(k) & D(B(m+1))>logStock(k)
                c(k)=polyval(Coef(m,:),logStock(k)-D(B(m)));
                flag=1;
                break;
            else
                ;
            end
        end
        %if flag==0
        if flag==0
            disp('not find the value for ')
            exp(logStock(k))*ScaleUnit   
        end
    end
end
val=ScaleUnit*c;
%t2=toc;
%T=T+t2;
%end
%T/10000=0.17623696944723

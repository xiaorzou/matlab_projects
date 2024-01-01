%function void = main(model, strike, fftpower,  skip1,skip2, deg,m, left, right)


function void = drive()
Model=1;
%strike=[90,95,100,105];
strike = [100];
fftpower=10;
deg=2;
m=5;
%alpha=[1.1,1.2,1.3,1.4,1.5,2];

left=-1;
right=1;
%numA=32*2;
%numB=32*2;
alpha=1;
beta=1;
numPeriods=10;
%results=zeros(6, 2);
%sprintf('alpha=%f,beta=%f,left=%f,right=%f,deg=%f,fftpower=%f',alpha,beta,left,right,deg,fftpower)
%if Model~=0
%model=Model;
for model=Model:Model
%for model=1:6
for s=1:length(strike)
    results=zeros(1,2);
    %fftpower=15;
    %alpha=1.2;
    numA=2^6;
    numB=2^6;
    for i=1:1
    %for fftpower=10:14
    %if model==1 
    %    left=-1.2;
    %    right=1.2;
    %elseif model==2
    %    left=-2;
    %    right=2;
    %elseif model==3 | model==4
    %%
    % 
    %   for x = 1:10
    %       disp(x)
    %   end
    % 
    %    left=-2.2;
    %    right=2.2;
    %else
    %    left=-2;
    %    right=2;
    %end
        tic;
        val=main(model, strike(s), fftpower, deg,m, left, right, numPeriods, alpha, beta, numA, numB);
        t=toc;
        results(i, 1)=val;
        results(i,2)=t;
        %alpha=alpha+0.2;
        fftpower=fftpower+1;
        numA=numA*2;
        numB=numB*2;
        %[numA, numB, fftpower]
    end
    sprintf('the model is %f.',model)
    results(:,1)
    results(:,2)
end
%results(5, model)=t;
%    sprintf('the model is %f.',model)
%results(:, model)
end

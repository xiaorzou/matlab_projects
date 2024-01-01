function Fib = Fib( ndim )
%{
if ndim= 4, output is
     1     1     0     0
     1     2     1     0
     1     3     3     1
%}
Fib=[1,0;1,1];
for i=3:(ndim)
    [rows, cols]=size(Fib);
    z1=zeros(rows,1);
    Fib=[Fib,z1];
    last=Fib(end,:);
    nr=zeros(1, cols+1);
    nr(1)=1;
    for i=1:cols
        nr(1+i)=last(i)+last(i+1);
    end
    Fib=[Fib;nr];
end
end


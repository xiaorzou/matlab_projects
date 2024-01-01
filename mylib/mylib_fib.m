function fib = mylib_fib( ndim )
    %{
    if ndim= 4, output is
         1     1     0     0
         1     2     1     0
         1     3     3     1
    %}
    fib=[1,0;1,1];
    for i=3:(ndim)
        [rows, cols]=size(fib);
        z1=zeros(rows,1);
        fib=[fib,z1];
        last=fib(end,:);
        nr=zeros(1, cols+1);
        nr(1)=1;
        for j=1:cols
            nr(1+j)=last(j)+last(j+1);
        end
        fib=[fib;nr];
    end
end


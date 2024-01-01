%{
get w^{(n)}(x): the n-th dervative values of w(x), n is in [0,5]
where w(x)=const*exp(k/((x+a)(x+b)) for x in (-b,-a) and (a,b) are in the
standard setting for fourier_normalizer project.
%}
function val = fourier_normalizer_get_w_der(x,a,b,c,d,n)
    if b < a
            print('b>a is required');
            exit(0);
    end
    if c < 0 || a < 0
            print('a, c>0 is required');
            exit(0);
    end

    if d <= 0
            print('d>0 is required');
            exit(0);
    end

    part1 = x(x<=-b);
    part2 = x(x>-b & x<-a);
    part3 = x(x>=-a & x<=a);
    part4 = x(x>a & x<b);
    part5 = x(x>= b);
    delta = b-a;
    k = delta^2*d/(4*c^2);
    val1 = zeros(1,length(part1));
    val2 = fourier_normalizer_get_w_der_leftwing(part2, a, b, n, k);
    val3 = zeros(1,length(part3));
    val4 = -fourier_normalizer_get_w_der_leftwing(-part4, a, b, n, k)*(-1)^n;
    val5 = zeros(1,length(part5));
    val = [val1,val2,val3,val4,val5];
end


%{
get w^{(n)}(x): the n-th dervative values of w(x) over (-b,-a),
where w(x)=const*exp(k/((x+a)(x+b)) for x in (-b,-a) and (a,b) are in the
standard setting for fourier_normalizer project.

input:
    @n: degree of derivative [0,5],  0 refers w(x) 
    @a,b,k: define the function exp(k/((x+a)(x+b)) used in w(x)
    @x: a point in (-b,-a)
output:
    @val: w^{(n)}(x)
%}
function val = fourier_normalizer_get_w_der_leftwing(x, a, b, n, k)
    y0 = get_der_special(x,a,b,0,k);
    if n == 0
        val = exp(y0);
    else
        y1 = get_der_special(x,a,b,1,k);
        if n==1
            val = exp(y0).*y1;
        else 
            y2 = get_der_special(x,a,b,2,k);
            if n ==2
                val = exp(y0).*(y1.^2+y2);
            else
                y3 = get_der_special(x,a,b,3,k);
                if n ==3
                    val = exp(y0).*(y1.^3+3*y1.*y2+y3);
                else
                    y4 = get_der_special(x,a,b,4,k);
                    if n==4
                        val = exp(y0).*(y1.^4+6*y1.^2.*y2+3*y2.^2+4*y1.*y3+y4);
                    else
                        y5 = get_der_special(x,a,b,5,k);
                        if n == 5
                            val = exp(y0).*(y1.^5+10*y1.^3.*y2+10*y1.^2.*y3+15*y1.*y2.^2+5*y1.*y4+10*y2.*y3+y5);
                        else
                            print(['not implmented.  n=', int2str(n)])
                        end
                    end  
                end
            end
        end        
    end
end


% f(x) = k/((b+x)(a+x))=k/delta(1/(x+a)-1/(x+b))
%output f^(n)(x)
function val = get_der_special(x,a,b,n,k)

    delta = (b-a);
    if n==0
        val = k*((x+a).^(-1)-(x+b).^(-1))/delta;
    else
        val = k*factorial(n)*((x+a).^(-(n+1))-(x+b).^(-(n+1)))*(-1)^n/delta;
    end
end


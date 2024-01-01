function val = normalizer(x,alpha,beta,a)
% u_ext(-x) = -u_ext(-x)  odd function
% u_ext(x) = 0 if |x|>=beta, or [-alpha, alpha]
% u_ext(x) = u(x) for x in (alpha, beta)
if beta < alpha
        print('beta>alpha is required');
        exit(0);
end
if a < 0 || alpha< 0
        print('a, alpha>0 is required');
        exit(0);
end

part1 = x(x<=-beta);
part2 = x(x>-beta & x<-alpha);
part3 = x(x>=-alpha & x<=alpha);
part4 = x(x>alpha & x<beta);
part5 = x(x>= beta);

val1 = zeros(1,length(part1));
val2 = -normalizer_core(-part2, alpha, beta, a);
val3 = zeros(1,length(part3));
val4 = normalizer_core(part4, alpha, beta, a);
val5 = zeros(1,length(part5));

val = [val1,val2,val3,val4,val5];
end




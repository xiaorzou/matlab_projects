function  val = normalizer_core(x,alpha,beta,a)
% x in (alpha, beta)
% u(x)=(beta-alpha)^3((alpha+beta)/2-x)exp((beta-alpha)^2/(4a^2(alpha-x)(beta-x)))/(2a^2(alpha-x)^2(beta-x)^2),
delta = beta-alpha;
expression = (alpha-x).*(beta-x);
val = -2*(a/delta)*(exp(delta^2./(4*a^2*expression)));
end


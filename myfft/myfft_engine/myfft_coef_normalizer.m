
%input:
%   @coef:  the coefficents of cos expansion of some density function q(x)
%   @h: normalier, it is the coefficents of cos expansition of normalizer,
%   which is equal to 1 over the range (-a,a), and 0 outside (-b,b). 
%   @n: convolution step over coef(1:n), i.e. n=dim(coef_normalized)
%output:
%   @coef_normalized:  the output array,  for 1<=k<=n
%       coef_normalized
%           
%     
%   
function coef_normalized = myfft_coef_normalizer(coef, h, n)
    M = length(h);
    N = length(coef);
    if n>(N-M)
        fprintf('%s\n','invalide n, n>length(coef)-length(h) is required');
        exit(0)
    end
    coef_normalized = zeros(1,n);
    %coef_fli = fliplr(coef(1:M));
    h_fli = fliplr(h);
    for k = 1:n
        I = dot(h,coef(k:k+M-1));
        if k<M
            %II = dot(h(1:k),coef_fli(M-k+1:M)) + dot(h(k+1:M), coef(2:M-k+1));
            II = dot(coef(1:k),h_fli(M-k+1:M)) + dot(coef(2:M-k+1),h(k+1:M));
            %II-II_new
        else
            II = dot(h_fli(1:M),coef(k-M+1:k));
        end
        coef_normalized(1,k)=(I+II)/2;
    end
    
    %add the followin two lines since they are removed in myfft_option_test_enhancement_convergence
    coef_normalized(1) = coef(1)*h(1) + 0.5*dot(h(2:M), coef(2:M)); %newly added, 
    coef_normalized(2:M) = coef_normalized(2:M)+ 0.5*coef(1).*h(2:M); %new added,  see myfft.pdf
end


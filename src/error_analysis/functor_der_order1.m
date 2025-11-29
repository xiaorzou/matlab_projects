function fun =functor_der_order1(model, para, b)
     if strcmp(model, 'power')
        fun = @(x) -(2/b^2)*para*(1-(x/b).^2).^(para-1).*x; %para=deg
    elseif strcmp(model, 'cos_alpha')
        fun = @(x) -(para*pi/b)*sin(para*x*(pi/b)); %para=alpha
    elseif strcmp(model, 'sin_alpha')
        fun = @(x) (para*pi/b)*cos(para*x*(pi/b)); %para=alpha
    elseif strcmp(model, 'ln')
        fun = @(x) (2*para/(log(b^2+1))^para)*(log(x^2+1)).^(para-1).*x./(x^2+1); %para=deg
    elseif strcmp(model, 'abs')
        fun = @(x) -para*(b-abs(x)).^(para-1).*(sign(x)+1-sign(x).^2); %para=deg
    elseif strcmp(model, 'exp')
        fun = @(x) (-para*(abs(x)).^(para-1)).*(sign(x)+1-sign(x).^2).*exp((b^para-abs(x).^para))/exp(b^para); %para=deg
    elseif strcmp(model, 'cos')
        fun = @(x)dot(c_true.*(0:(para-1)), -sin(((0:(para-1))*x))); %para= terms
    else
        fprintf('invalid model %s', model)
    end   
end


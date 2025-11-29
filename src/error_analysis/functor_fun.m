function fun =functor_fun(model, para, b)
    if strcmp(model, 'power')
        fun = @(x)(1-(x/b).^2).^para; %para=deg
    elseif strcmp(model, 'cos_alpha')
        fun = @(x)cos((para*pi/b)*x); %para=alpha
    elseif strcmp(model, 'sin_alpha')
        fun = @(x)sin((para*pi/b)*x); %para=alpha
    elseif strcmp(model, 'ln')
        fun = @(x)(log(x.^2+1)).^para/(log(b^2+1))^para; %para=deg
    elseif strcmp(model, 'abs')
        fun = @(x)(b-abs(x)).^para; %para=deg
    elseif strcmp(model, 'exp')
        fun = @(x)exp((b^para-abs(x).^para))/exp(b^para); %para=deg
    elseif strcmp(model, 'cos')
        fun = @(x)dot(c_true, cos(((0:(para-1)).*x))); %para= terms
    else
        fprintf('invalid model %s', model)
    end
end

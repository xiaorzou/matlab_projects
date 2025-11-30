
function fun =functor_der_order2(model, para, b)
     if strcmp(model, 'power')
        if para == 1
            fun = @(x) -(2/b^2)*ones(1,length(x)); %para=deg
        else
            fun = @(x) (4/b^4)*para*(para-1).*(1-(x/b).^2).^(para-2).*x.^2  -(2/b^2)*para*(1-(x/b).^2).^(para-1) ; %para=deg
        end
    elseif strcmp(model, 'cos_alpha')
        fun = @(x) -(para*pi/b)^2*cos(para*x); %para=alpha
    elseif strcmp(model, 'sin_alpha')
        fun = @(x) -(para*pi/b)^2*sin(para*x); %para=alpha
    elseif strcmp(model, 'ln')
        if para==1
            fun = @(x) (2*para/(log(b^2+1))^para).*(1-x.^2)./(x^2+1).^2; %para=deg
        else
            fun = @(x) (4*para*(para-1)/(log(b^2+1))^para)*(log(x^2+1)).^(para-1)*x^2./(x^2+1).^2 + (2*para/(log(b^2+1))^para)*(log(x^2+1)).^(para-1).*(1-x.^2)./(x^2+1).^2; 
        end
    elseif strcmp(model, 'abs')
        if para == 1
            fun = @(x) zeros(1,length(x));
        else
            fun = @(x) para*(para-1)*(b-abs(x)).^(para-2).*(sign(x)+1-sign(x).^2); %para=deg
        end
    elseif strcmp(model, 'exp')
        if para == 1
            fun = @(x) (sign(x)+1-sign(x).^2).*exp((b^para-abs(x)))/exp(b^para); %para=deg
        else
            fun = @(x) (-para*(para-1)*(abs(x)).^(para-2)).*(sign(x)+1-sign(x).^2).*exp((b^para-abs(x).^para))/exp(b^para) + (para^2*(x.^2).^(para-1)).*exp((b^para-abs(x).^para))/exp(b^para); %para=deg
        end
    elseif strcmp(model, 'cos')
        fun = @(x)dot(c_true.*(0:(para-1)).^2, -cos(((0:(para-1))*x))); %para= terms
    else
        fprintf('invalid model %s', model)
    end   
end



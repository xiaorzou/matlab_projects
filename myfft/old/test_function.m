function val = test_function(x, opt)
if opt==1
    val = x;
elseif opt == 2
    val = x.^3;
elseif opt == 3
    val = exp(x) - exp(-x);
elseif opt == 4
    val = sin(x)-sin(2*x)+sin(10*x);
end


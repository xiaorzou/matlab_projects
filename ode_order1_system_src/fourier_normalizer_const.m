
function fourier_normalizer_const()
    global cuf_off_para normalizer_d normalizer_q options
    global test_type weight_constrain s e dim ps qs c3 degs init_method qqs  paras ws pps
    global maxfunevals figtype_fig plot_location flag_plot flag_applying_constrains
    flag_display = true;
    cuf_off_para = 0.5;
    normalizer_d = 20; %default value of parameter d for h 
    normalizer_q = 10; %default value of parameter q for h
    maxfunevals = 1000;%setting for paper
    maxiter = maxfunevals; %setting for paper
    constrainttolerance = 10^(-20); %setting for paper
    constrainttolerance_tolx = 10^-15; %setting for paper

    flag_plot = false;
    not_default = true;
    figtype_fig = '.fig';
    test_type = 'base';  
    plot_location = 'northwest';    

    s = 1; %do not change
    e = 3*s; %do not change
    dim = 3; %do not change  
    
    weight_constrain = 1;%setting in doc
    ps = [0.1,0.1,0.1]; %setting in doc
    qs = [0.1,0.1,0.1]; %setting in doc
    c3 = 1;
    flag_applying_constrains = [true, false]; %setting in doc
    degs = [0,1,2]; %setting in doc
    paras = [1,5]; %setting in doc
    qqs = [6,7,8]; %setting in doc
    init_method = 'rk'; %setting in doc
    
    if not_default %define specific test in following section!
        constrainttolerance = 10^(-25); %setting for paper
        constrainttolerance_tolx = 10^-25; %setting for paper
        maxfunevals = 10000;
        %constrainttolerance_tolx = 10^-20; 
        %ps = [0.0,0.0,0.0]; %setting in complete testing
        flag_applying_constrains = [true]; 
        qqs = [6];
        paras = [1,5];
        degs = [2];    
        %maxfunevals = 10000;
        %init_method = 'simple';
        %init_method = 'cheat';
    end
    
    %init_method = 'cheat';
    %init_method = 'const';
    %init_method = 'simple'; %use linear method to guess init
    
    ws = ones(1,dim); %weights for optimization
    pps = ones(1,dim); % p values of coefficients in dynamics,  so we have the flexibity to handle algebra equation

    if flag_display
        options = optimoptions(@fmincon, ...
            'Algorithm', 'interior-point',...
            'Display', 'Iter-detailed',...
            'GradObj','on',...
            'MaxIter', maxiter,...
            'TolFun', constrainttolerance,...
            'TolX', constrainttolerance_tolx,...
            'TolCon', constrainttolerance,...
            'ObjectiveLimit',-1e30, ...
            'MaxFunEvals', maxfunevals,...
            'GradConstr','on');
    else
        options = optimoptions(@fmincon, ...
            'Algorithm', 'interior-point',...
            'GradObj','on',...
            'MaxIter', maxiter,...
            'TolFun', constrainttolerance,...
            'TolX', constrainttolerance_tolx,...
            'TolCon', constrainttolerance,...
            'ObjectiveLimit',-1e30, ...
            'MaxFunEvals', maxfunevals,...
            'GradConstr','on');
    end
   
end
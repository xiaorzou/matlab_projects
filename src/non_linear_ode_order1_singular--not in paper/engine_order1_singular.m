
%{
it solves general ode y'(x) = fun(x,y) 
we solve u'(x)= fun(x,u)*h(x) with h(x) as cut-off function
 inputs   
    fun: function para for right of ODE
    partial_u: function para for f'_u(x,u)*h(x), function parameter for partial dervative to y
    
important note:  both fun and partial_u should be defined in absolute
sense.   i.e. they are evenly extended to left.  fun(x,y)=f(-x,y),
partial_u(-x,y)= partial_u(x,y).  we manually time -1 to deal with z
values in the implemetation.

    y_s:  = y(s)=y(-s)
    z_s:  = y'(s)=-y'(-s),  which should be able to get from ODE itself. need not be accurate, we use it to do initial guess. 
    s,e,p,q:  four parameters in construction for h(x)
output:
    x_N:  x values at grid points
    u_N:  u values at x_N
    z_N:  u'(x) values at x_N
    pos:  index positions for grid points in [s,e]
%}
function [z_N, u_N, myobjopt] = engine_order1_singular(fun, partial_u, y_s, z_s, s, e, p,q, G_L, R, eta, z_init, exclude_x0)
    %addpath('D:/matlab/mylib') 
    
    
    flag_display = false;
    maxiter = 1000;
    maxfunevals = maxiter;
    constrainttolerance = 10^-10;
    
    fourier_normalizer_const()
    %s = 1;
    %e = 3*s;
    %q = 9;
    %p = q-1; %p should be smaller than q,  
    M = 2^q;
    N = 2*M;
    n = 2^p; % we should make n as large as possible,  as such, p=q-1 should be always true
    lambda = (e-s)/n;
    m = (M-n)/2;
    delta = lambda*m;
    o = s - delta;
    b = e-s + 2*delta;
    x_N = -b + o + lambda*(0:N-1);
    x_R = x_N(M+1:end);

    
    J_M = zeros(1,M);
    J_M(2:M) = 1./(1:M-1);
    
    J_N = zeros(1,N);
    J_N(1:M) = J_M;
    
    Phi_M = J_M.* cos(2*pi*(0:M-1)*(m+n)/N);
    Phi_N =  zeros(1,N);
    Phi_N(1:M) = Phi_M;
    bopi = b/pi;
    ifft_Phi = imag(ifft(Phi_N));
    

 

    
    
    
    u_N = zeros(1,N);
    %plot(x_M, z_init(1:M), x_M, z_true);
    %plot(x_M, u_true, x_M, u_init);
    
    zlast = [];
    myobjlast = [];
    glast = [];
    hlast = [];
    dmyobjlast = [];
    dglast = [];
    dhlast = [];
    lb = [];
    ub = [];
    
    function [myobj,g,h,dmyobj,dg,dh] = objcon(z) 
        %plot(x_M, z_M, x_M, z_true);
        if exclude_x0
            z_M =  [0, z(1:m+n-1),-z_s, z(m+n:end)];
        else
            z_M =  [z(1:m+n),-z_s, z(m+n+1:end)];
        end
        z_N = [z_M,0, -fliplr(z_M(2:M))];     
        a_0 = y_s + 2*bopi * dot(Phi_N, imag(ifft(z_N)));
        u_N = a_0 - 2*bopi*N*real(ifft(J_N.*(imag(ifft(z_N))))); % we need u_N, rather than just u_M, u(0) is not known!
        xu_R = [x_R;u_N(M+1:N)];
        F_R = fun(xu_R, R); 
        F_M = [0,-fliplr(F_R(2:M))];
        myobj = sum(eta.*(G_L.*z_M-F_M).^2)/N;
        myobj = myobj';
        DF_R = partial_u(xu_R);
        %DF_M = -partial_u(xu_R);
        DF_M = [0,-fliplr(DF_R(2:M))];
        Psi_M = eta.*(G_L.*z_M-F_M).*DF_M;
        I = sum(Psi_M);
        Psi_N = zeros(1,N);
        Psi_N(1:M) = Psi_M;
        W = 4*bopi*I*ifft_Phi - 4*bopi*N*imag(ifft(J_N.*real(ifft(Psi_N))));
        dmyobj = (eta.*(G_L.*z_M-F_M)-W(1:M))/M;
        dmyobj = dmyobj(1:M);
        if exclude_x0
            dmyobj = [dmyobj(2:m+n),dmyobj(m+n+2:end)]; 
        else
            dmyobj = [dmyobj(1:m+n),dmyobj(m+n+2:end)]; 
        end
        dmyobj = dmyobj';
        
        g = [];
        h = [];
        dg = [];
        dh = [];
    end
    if flag_display
        options = optimoptions(@fmincon, ...
            'Algorithm', 'interior-point',...
            'Display', 'Iter-detailed',...
            'GradObj','on',...
            'MaxIter', maxiter,...
            'TolFun', constrainttolerance,...
            'TolX', constrainttolerance,...
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
            'TolX', constrainttolerance,...
            'TolCon', constrainttolerance,...
            'ObjectiveLimit',-1e30, ...
            'MaxFunEvals', maxfunevals,...
            'GradConstr','on');
    end
        

    function [myobj,dmyobj] = obj(z)
        if ~isequal(z,zlast)
            [myobjlast, glast, hlast, dmyobjlast, dglast, dhlast] = objcon(z);
            zlast = z;
        end
        myobj = myobjlast;
        dmyobj = dmyobjlast;
    end
    function [g,h,dg,dh] = con(z)
        if ~isequal(z,zlast)
            [myobjlast, glast, hlast, dmyobjlast, dglast, dhlast] = objcon(z);
            zlast = z;
        end
        g = glast;
        h = hlast;
        dg = dglast;
        dh = dhlast; 
    end
    [zopt, myobjopt,iter, opt] = fmincon(@obj, z_init, [],[], [], [], lb, ub, @con, options);
    if exclude_x0
        zopt =  [0, zopt(1:m+n-1),-z_s, zopt(m+n:end)];
        %z_init = [0, z_init(1:m+n-1),-z_s, z_init(m+n:end)];  %make sense
        %if z_init is true value
    else
        zopt =  [zopt(1:m+n),-z_s, zopt(m+n+1:end)];
        %z_init = [z_init(1:m+n),-z_s, z_init(m+n+1:end)];
    end
    %max(abs(zopt(m+1:m+n+1)-z_init(m+1:m+n+1)))
    z_N = [zopt, 0, -fliplr(zopt(2:M))];
end

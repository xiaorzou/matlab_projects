
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
function [x_N, u_N , u_init_N, z_N, z_init_N, pos, myobjopt] = fourier_normalizer_app_ode_deg_1_general(fun, partial_u, y_s, z_s, s, e, p,q, init_method)
    addpath('D:/matlab/fourier_normalizer/fourier_normalizer_engine')  
    addpath('D:/matlab/cos_approx/cos_approx_engine')  
    %addpath('D:/matlab/mylib') 
    global cuf_off_para 
    
    flag_display = false;
    maxiter = 1000;
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
    x_M = x_N(1:M);
    x_R = x_N(M+1:end);
    co_R = fourier_normalizer_cut_off(x_R, s-delta,s,e,e+delta,cuf_off_para);
    co_M = [0,fliplr(co_R(2:M))];
    %co_M = [0,fliplr(co_R(2:M))];
    %co_even = co_M;
    %co_M = -co_M; %we only do calcuation on left!
    
    
    %{
    function val = fun(X)
        val = (abs(X(1,:)).^2 -abs(X(1,:)).^2.*X(2,:)).*co_M;
    end

    function val = fun_u(X)
        val = -abs(X(1,:)).^2.*co_M;
    end
    %}
  
    
    J_M = zeros(1,M);
    J_M(2:M) = 1./(1:M-1);
    
    J_N = zeros(1,N);
    J_N(1:M) = J_M;
    
    Phi_M = J_M.* cos(2*pi*(0:M-1)*(m+n)/N);
    Phi_N =  zeros(1,N);
    Phi_N(1:M) = Phi_M;
    bopi = b/pi;
    ifft_Phi = imag(ifft(Phi_N));
    

 
    if strcmp(init_method, 'simple')
       %u_init = u_true;
        z_init = zeros(1,M);
        u_init = zeros(1,M);
        index_0 = m+n+1;
        z_init(index_0) = z_s;  
        u_init(index_0) = y_s; 
        for i = (index_0+1):M
                u_init(i) = (z_init(i-1)*lambda + u_init(i-1))*co_M(i); %u is even!
                z_init(i) = -fun([x_M(i); u_init(i)])*co_M(i); %z is odd!
        end
        for i = (index_0-1):-1:2
                u_init(i) = (z_init(i+1)*(-lambda) + u_init(i+1))*co_M(i); %u is even!
                z_init(i) = -fun([x_M(i); u_init(i)])*co_M(i); %z is odd!
        end
    else
        u_rk = runge_kutta_ode_orde_1(fun, x_R,y_s, m+1, []);
        u_init =  [0,fliplr(u_rk(2:M))];
        u_init = u_init.*co_M;
        z_init =  -fun([x_M; u_init]).*co_M;
    end
    
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
    
    function [myobj,g,h,dmyobj,dg,dh] = objcon(z_M) 
        %plot(x_M, z_M, x_M, z_true);
        z_N = [z_M,0, -fliplr(z_M(2:M))];     
        a_0 = y_s + 2*bopi * dot(Phi_N, imag(ifft(z_N)));
        u_N = a_0 - 2*bopi*N*real(ifft(J_N.*(imag(ifft(z_N))))); % we need u_N, rather than just u_M, u(0) is not known!
        xu_M = [x_M;u_N(1:M)];
        F_M = -fun(xu_M).*co_M; 
        myobj = sum((z_M-F_M).^2)/N;
        myobj = myobj';
        DF_M = -partial_u(xu_M).*co_M;
        Psi_M = (z_M-F_M).*DF_M;
        I = sum(Psi_M);
        Psi_N = zeros(1,N);
        Psi_N(1:M) = Psi_M;
        W = 4*bopi*I*ifft_Phi - 4*bopi*N*imag(ifft(J_N.*real(ifft(Psi_N))));
        dmyobj = (z_M-F_M-W(1:M))/M;
        dmyobj = dmyobj(1:M);
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
            'GradConstr','on');
    else
        options = optimoptions(@fmincon, ...
            'Algorithm', 'interior-point',...
            'GradObj','on',...
            'MaxIter', maxiter,...
            'TolFun', constrainttolerance,...
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
    pos = find(x_N<=e & x_N>=s);
    %u_R = u_N(M+1:end);
    %z_R = [0, -fliplr(zopt(2:M))];
    z_N = [zopt, 0, -fliplr(zopt(2:M))];
    %u_init_R = [0, fliplr(u_init(2:M))];
    u_init_N = [u_init, 0, fliplr(u_init(2:M))];
    %z_init_R = [0, -fliplr(z_init(2:M))];
    z_init_N = [z_init, 0, -fliplr(z_init(2:M))];
    
end

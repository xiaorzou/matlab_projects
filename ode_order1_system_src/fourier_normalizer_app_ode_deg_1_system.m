
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

    y_s_vect:  values at s
    s,e,p,q=p-1:  four parameters in construction for h(x)
output:
    u_N_out:  u values at x_N
    z_N_out:  u'(x) values at x_N
    coef:  coefs for u
%}

function [u_N_out, z_N_out, coef, fval] = fourier_normalizer_app_ode_deg_1_system(fun, partial_u, y_s_vect, s, e, p, init_guess)
    addpath('D:/matlab/fourier_normalizer/fourier_normalizer_engine')  
    addpath('D:/matlab/cos_approx/cos_approx_engine')  
    flag_display = false;
    maxiter = 1000;
    constrainttolerance = 10^-10;
    fourier_normalizer_const()
    
    [M, dim] = size(init_guess);
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
    
    zlast = [];
    myobjlast = [];
    glast = [];
    hlast = [];
    dmyobjlast = [];
    dglast = [];
    dhlast = [];
    lb = [];
    ub = [];
    
    % z_comb is for the left part:  z_comb(1)..., z_comb(M) are assoicte to
    % the first half of z vector with N element), i.e. they are assoicated
    % to x_1=-b, ..., x_M=-lambda, note that x_{M+1}=0. 
    function [myobj_comb,g,h,dmyobj_comb,dg,dh] = objcon(z_comb)
        dmyobj_comb = zeros(dim*M,1);
        myobj_comb = 0;
        xu_comb_right = zeros(M, dim+1); % for the righ half !!!
        xu_comb_right(:,1) = x_R;
        for alpha = 1:dim
            z_M = z_comb((M*(alpha-1)+1):M*alpha);
            z_N = [z_M,0, -fliplr(z_M(2:M))];     
            a_0 = y_s_vect(alpha) + 2*bopi * dot(Phi_N, imag(ifft(z_N)));
            u_N = a_0 - 2*bopi*N*real(ifft(J_N.*(imag(ifft(z_N))))); % we need u_N, rather than just u_M, u(0) is not known!
            xu_comb_right(:,alpha+1) = u_N(M+1:end)';
        end
        
        F_M_comb_right = fun(xu_comb_right);  %the values on the right half!!!
        
        for alpha = 1:dim %index of u, i.e. alpha in doc
            Psi_M_alpha_left = zeros(M,1);
            for beta = 1:dim %index of F, i.e. beta in doc
                DF_M_right = partial_u(xu_comb_right, beta, alpha);
                DF_M_left = -[0;fliplr(DF_M_right(2:M)')'];
                z_M_left = z_comb((M*(beta-1)+1):M*beta);
                F_M_right = F_M_comb_right(:,beta);
                F_M_left = -[0;fliplr(F_M_right(2:M)')'];
                Psi_M_alpha_left = Psi_M_alpha_left + (z_M_left'-F_M_left).*DF_M_left;
            end
            I =  sum(Psi_M_alpha_left);
            Psi_N_alpha = zeros(1,N);
            Psi_N_alpha(1:M) = Psi_M_alpha_left;
            F_M_right = F_M_comb_right(:,alpha);
            F_M_left = -[0;fliplr(F_M_right(2:M)')'];
            z_M_left = z_comb((M*(alpha-1)+1):M*alpha);
            myobj_comb = myobj_comb + sum((z_M_left-F_M_left').^2)/(dim*N);
            
            W = 4*bopi*I*ifft_Phi - 4*bopi*N*imag(ifft(J_N.*real(ifft(Psi_N_alpha))));
            dmyobj_comb((M*(alpha-1)+1):M*alpha)= (z_M_left-F_M_left'-W(1:M))/(dim*M);
        end
        %myobj_comb = myobj_comb';
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
    z_init = zeros(1, dim*M);
    for alpha_ = 1:dim
        z_init((M*(alpha_-1)+1):M*alpha_) = init_guess(:,alpha_);
    end
    [zopt, fval] = fmincon(@obj, z_init, [],[], [], [], lb, ub, @con, options);
    %pos = find(x_N<=e & x_N>=s);
    %u_R = u_N(M+1:end);
    %z_R = [0, -fliplr(zopt(2:M))];
    z_N_out = zeros(N, dim);
    u_N_out = zeros(N, dim);
    coef = zeros(M,dim);
    for a = 1:dim
        z = zopt((M*(a-1)+1):M*a);
        z = [z,0, -fliplr(z(2:M))];     
        a_0_this = y_s_vect(a) + 2*bopi * dot(Phi_N, imag(ifft(z)));
        u = a_0_this - 2*bopi*N*real(ifft(J_N.*(imag(ifft(z))))); 
        z_N_out(:, a)= z;
        u_N_out(:, a) = u;
        coef(:, a) = cos_approx_engine_value2coef_wo_loop(u);
    end
end

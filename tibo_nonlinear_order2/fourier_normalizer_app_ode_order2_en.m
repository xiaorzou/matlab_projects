
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
%fourier_normalizer_app_ode_deg_2_general(@fun, @partial_u,@partial_v, y_s, z_s, s, e, p,q, init_method, alpha, beta, cond_type)

%function [x_N, z_N, b,o,m,n, a_0, a_1, v_rk,u_rk,z_rk, pos, myobjopt] = fourier_normalizer_app_ode_order2_en(fun, partial_u, partial_v, s, e, p,q, alpha, beta, cond_type, tor_bvp, counter_bvp, use_bvp, init_slope)
%function [X_inte, X_rk4, a_0, a_1, b,o, myobjopt_inte, myobjopt_shooting] = fourier_normalizer_app_ode_order2_en(fun, partial_u, partial_v, s, e, p,q, greeks, AA, init_guess, flag_display,...
function [X_inte, a_0, a_1, b,o, myobjopt_inte] = fourier_normalizer_app_ode_order2_en(fun, partial_u, partial_v, s, e, p,q, greeks, AA, flag_display,...
    constrainttolerance,maxiter, apply_condition_diri,cond_b, ...
    apply_condition_on_lower_boundary_flag, apply_condition_on_lower_boundary_value, apply_condition_on_upper_boundary_flag, apply_condition_on_upper_boundary_value, z_init)
    addpath('D:/matlab/fourier_normalizer/fourier_normalizer_engine')  
    addpath('D:/matlab/cos_approx/cos_approx_engine')  
    %addpath('D:/matlab/mylib') 
    global cuf_off_para 
    
    %flag_display = false;
    %flag_display = true;
    %maxiter = 10000;
    %constrainttolerance = 10^-10;
    
    
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
    co_R = fourier_normalizer_cut_off(x_R, s-delta,s,e,e+delta,cuf_off_para);
    
    function val = fun_grid(X)   
        val = fun(X);
        if length(val)==length(co_R)
            val = val.*co_R;
        else
            x1 = X(1,:);
            t=length(val);
            for i=1:t
                pos_this = find(x_R>=x1(i));
                val(i) =val(i)*co_R(pos_this(1));
            end
        end
    end
    function val = partial_u_grid(X)
        val = partial_u(X);
        if length(val)==length(co_R)
            val = val.*co_R;
        else
            t=length(val);
            x1 = X(1,:);
            for i=1:t
                pos_this = find(x_R>=x1(i));
                val(i) =val(i)*co_R(pos_this(1));
            end
        end
    end
    function val = partial_v_grid(X)
        val = partial_v(X);
        if length(val)==length(co_R)
            val = val.*co_R;
        else
            t=length(val);
            x1 = X(1,:);
            for i=1:t
                pos_this = find(x_R>=x1(i));
                val(i) =val(i)*co_R(pos_this(1));
            end
        end
    end

    J_M = zeros(1,M);
    J_M(2:M) = 1./(1:M-1);
    
    
    
    J_N = zeros(1,N);
    J_N(1:M) = J_M;
    A = ones(1,N);
    A(2:2:end) = -1;
    
      
    Phi_M_cmpn = cos(2*pi*(0:M-1)*(m+n)/N);
    Phi_N_cmpn =  zeros(1,N);
    Phi_N_cmpn(1:M) = Phi_M_cmpn;
    Phi_M_smpn = sin(2*pi*(0:M-1)*(m+n)/N);
    Phi_N_smpn =  zeros(1,N);
    Phi_N_smpn(1:M) = Phi_M_smpn;
    
    Phi_M_s = sin(2*pi*(0:M-1)*(m)/N);
    Phi_M_c = cos(2*pi*(0:M-1)*(m)/N);
    Phi_N_s =  zeros(1,N);
    Phi_N_s(1:M) = Phi_M_s;
    Phi_N_c =  zeros(1,N);
    Phi_N_c(1:M) = Phi_M_c;   
    bopi = b/pi;
    A_opt =[];
    b_opt = [];
    if all(AA==[[1,0,0,0];[0,0,1,0]]) 
        if strcmp(apply_condition_diri,'upbound') || strcmp(apply_condition_diri,'lowbound')  %'Dirichlet' only! %'Dirichlet' only!
            A2d = zeros(M,M);
            for ii = 1:M
                A2d(ii,:) = sin(2*pi*(M:N-1)*(ii-1)/N);
            end
            A_opt = (J_M.*Phi_M_c+J_M.^2.*Phi_M_s*bopi/(e-s) - J_M.^2.*Phi_M_smpn*bopi/(e-s));
            alt_Mm1 = ones(1,M);
            alt_Mm1(2:2:end)=-1;
            A_opt =A_opt.*alt_Mm1; 
            A_opt = (4/N)*(A_opt*A2d);
            if strcmp(apply_condition_diri,'upbound')
                A_opt = -A_opt;
                %b_opt = -(greeks(2)-greeks(1))/((e-s)*bopi);
                b_opt = -((greeks(2)-greeks(1))/(e-s)-cond_b)/bopi;
            else
                b_opt =((greeks(2)-greeks(1))/(e-s)-cond_b)/bopi;
            end
        end
    end
    Aeq_opt = [];
    beq_opt = [];
    ub_opt = [];
    lb_opt = [];
    nonlcon_opt = [];
    options_opt = [];
    

   %[init_opt, myobjopt_shooting, X_rk4] = fourier_normalizer_app_ode_get_init_value_by_rk4_order2(fun, s, e,  greeks, AA, n, init_guess, delta, maxiter, constrainttolerance);    
       
   %alpha =init_opt(1);
   %beta = init_opt(2);
   %if strcmp(cond_type, 'Neumann')
   %alpha_n =init_opt(2);
   %beta_n = init_opt(1);
   %else
   %alpha_d =init_opt(1);
   %beta_d = init_opt(2);
   %end

    m_11 = (AA(1,1)*s+AA(1,2)+AA(1,3)*e+AA(1,4));
    m_12 = (AA(1,1)+AA(1,3));
    m_21 = (AA(2,1)*s+AA(2,2)+AA(2,3)*e+AA(2,4));
    m_22 = (AA(2,1)+AA(2,3));
    delta = m_11*m_22 - m_12*m_21;
    grade_sm = 4*bopi^2*imag(ifft(A.*(J_N.^2).*Phi_N_s));
    grade_cm = 4*bopi*imag(ifft(A.*(J_N).*Phi_N_c));
    grade_smpn = 4*bopi^2*imag(ifft(A.*(J_N.^2).*Phi_N_smpn));
    grade_cmpn = 4*bopi*imag(ifft(A.*(J_N).*Phi_N_cmpn));
    %m_22*grad_I - m12*grad_II
    grad_1 = AA(1,1)*grade_sm + AA(1,2)*grade_cm+AA(1,3)*grade_smpn+AA(1,4)*grade_cmpn;
    grad_2 = AA(2,1)*grade_sm + AA(2,2)*grade_cm+AA(2,3)*grade_smpn+AA(2,4)*grade_cmpn;
    a_0_grad = 1/delta*(m_22*grad_1 - m_12*grad_2);
    a_1_grad = 1/delta*(-m_21*grad_1 + m_11*grad_2);
    a_0_grad = a_0_grad(M+1:N);
    a_1_grad = a_1_grad(M+1:N);
    %grade for Nueman 
    %{
    a_0_grad_old = 4*bopi*imag(ifft(A.*J_N.*Phi_N_c));
    a_0_grad_old = a_0_grad_old(M+1:N);
    a_1_grad_old = 4*bopi^2*imag(ifft(A.*(J_N.^2).*Phi_N_s));
    a_1_grad_old = a_1_grad_old(M+1:N);
    a_1_grad_old = a_1_grad_old - s*a_0_grad_old;
    %}
    %close grade for Dirichlet
    
    a_0_grad_old = 4*bopi^2/(e-s)*imag(ifft(A.*(J_N.^2).*(Phi_N_smpn-Phi_N_s)));
    a_0_grad_old = a_0_grad_old(M+1:N);
    a_1_grad_old = 4*bopi^2*imag(ifft(A.*(J_N.^2).*Phi_N_s));
    a_1_grad_old = a_1_grad_old(M+1:N);
    a_1_grad_old = a_1_grad_old - s*a_0_grad_old;
    
 
    %[~, z_rk] = runge_kutta_ode_order_k(@fun_grid, x_R,init_opt, m+1, []); 
    %z_init = z_rk;
    %if ~isempty(X_inte_n)
        %z_init = X_inte_n(4,:);
    %end
    
    %[X, z] = runge_kutta_ode_order_2(@fun_grid, x_R,inits, m+1, []); 
  
   

    zlast = [];
    myobjlast = [];
    glast = [];
    hlast = [];
    dmyobjlast = [];
    dglast = [];
    dhlast = [];
    
    
    lb = [];
    if apply_condition_on_lower_boundary_flag
        lb = apply_condition_on_lower_boundary_value *ones(M,1);
    end
    ub = [];
    if apply_condition_on_upper_boundary_flag;
      ub = apply_condition_on_upper_boundary_value*ones(M,1);  
    end
    a_0 = 0;
    a_1 = 0;
    u_R = [];
    v_R = [];
   if flag_display
        options = optimoptions(@fmincon, ...
            'Algorithm', 'interior-point',...
            'Display', 'Iter-detailed',...
            'GradObj','on',...
            'MaxIter', maxiter,...
            'MaxFunEvals',maxiter,...
            'TolFun', constrainttolerance,...
            'GradConstr','on');
    else
        options = optimoptions(@fmincon, ...
            'Algorithm', 'interior-point',...
            'GradObj','on',...
            'MaxIter', maxiter,...
            'MaxFunEvals',maxiter,...
            'TolFun', constrainttolerance,...
            'GradConstr','on');
    end
    %counter = 1;
    %record = zeros(maxiter,1);
    
    function [myobj,g,h,dmyobj,dg,dh] = objcon(z_M) 
        z_N = [0,-fliplr(z_M(2:M)), z_M];  
        
        ifft_this = imag(ifft(z_N));
        
        %if strcmp(cond_type, 'Neumann')
            %a_0_ne = alpha_n + 2*bopi*(sum(A.*J_N.*Phi_N_c.*ifft_this));
            %a_1_ne = beta_n -a_0_ne*s + 2*bopi^2*(sum(A.*(J_N.^2).*Phi_N_s.*ifft_this));
        %else
        %    %dict
            %a_0_d = (beta_d-alpha_d)/(e-s) + 2*bopi^2/(e-s)*(sum(A.*(J_N.^2).*(Phi_N_smpn-Phi_N_s).*ifft_this));
            %a_1_d = alpha_d -a_0_d*s + 2*bopi^2*(sum(A.*(J_N.^2).*Phi_N_s.*ifft_this));
        %end
        
        s_m = 2*bopi^2*(sum(A.*(J_N.^2).*Phi_N_s.*ifft_this));
        s_mn = 2*bopi^2*(sum(A.*(J_N.^2).*Phi_N_smpn.*ifft_this));
        c_m = 2*bopi*(sum(A.*(J_N).*Phi_N_c.*ifft_this));
        c_mn = 2*bopi*(sum(A.*(J_N).*Phi_N_cmpn.*ifft_this));

        b_1 = greeks(1)+ AA(1,1)*s_m + AA(1,2)*c_m + AA(1,3)*s_mn + AA(1,4)*c_mn;
        b_2 = greeks(2)+ AA(2,1)*s_m + AA(2,2)*c_m + AA(2,3)*s_mn+ AA(2,4)*c_mn;
        temp = [[m_11,m_12];[m_21,m_22]]\[b_1;b_2];
        a_0 = temp(1);
        a_1 = temp(2);
        
        
        u_R = a_0 - 2*bopi*N*real(ifft(J_N.*ifft_this)); % we need u_N, rather than just u_M, u(0) is not known!
        u_R = u_R(M+1:N);
        v_R = - 2*bopi^2*N*imag(ifft(J_N.^2.*ifft_this));
        v_R = v_R(M+1:N);
        v_R = a_1 + a_0* x_R +v_R;
        
        xvu_M = [x_R; v_R; u_R];
        F_M = fun_grid(xvu_M); 
        myobj = sum((z_M-F_M).^2)/N;
        myobj = myobj';
        %record(counter)=myobj;
        %counter = counter +1;
        %handle u=y'
        DF_u_M = partial_u_grid(xvu_M);
        Psi_u_M = (z_M-F_M).*DF_u_M;
        Psi_u_N = zeros(1,N);
        Psi_u_N(M+1:N) = Psi_u_M;
        I_u = sum (Psi_u_M);
        phi_u = 4*N*bopi*imag(ifft(J_N.*real(ifft(Psi_u_N))));
        phi_u = I_u*a_0_grad - phi_u(M+1:N);
        
        %handle v=y
        DF_v_M = partial_v_grid(xvu_M);
        Psi_v_M = (z_M-F_M).*DF_v_M;
        Psi_v_N = zeros(1,N);
        Psi_v_N(M+1:N) = Psi_v_M;
        I_v = sum(Psi_v_M);
        II_v = sum(x_R.*Psi_v_M);
        phi_v = 4*N*bopi^2*imag(ifft(J_N.^2.*imag(ifft(Psi_v_N))));
        phi_v = I_v*a_1_grad + II_v*a_0_grad -phi_v(M+1:N);
        
        %get gradient 
        W = phi_u + phi_v;
        dmyobj = (z_M-F_M-W)/M;
        dmyobj = dmyobj(1:M);
        dmyobj = dmyobj';
        g = [];
        h = [];
        dg = [];
        dh = [];
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
    %[init_opt, myobjopt_inte,~, ~] = fmincon(@obj, z_init, [],[], [], [], lb, ub, @con, options);
    [opt, myobjopt_inte,~, ~] = fmincon(@obj, z_init, A_opt,b_opt, [], [], lb, ub, @con, options);
    X_inte = [x_R; v_R; u_R; opt];
    %X_inte = [0,-fliplr(init_opt(2:M)), opt];
    %record;
end

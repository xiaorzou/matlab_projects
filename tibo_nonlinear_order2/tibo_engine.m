
%{

%}
function [X_inte, a_0, a_1, b,o, myobjopt_inte] = tibo_engine(fun, partial_u, partial_v, ...
    s, e, p,q, greeks, AA, z_init,r_values, ...
    condi_on_grid, condi_Diri, condi_v_low, condi_v_up, condi_u_low, condi_u_up, condi_z_low,  condi_z_up, up_val, low_val)
    %{
    condi_Diri,cond_b,...
    condi_v_low, condi_v_low_val, condi_v_up, condi_v_up_val, ...
    condi_u_low, condi_u_low_val, condi_u_up, condi_u_up_val, ...
    condi_z_low, condi_z_low_val, condi_z_up, condi_z_up_val)
    %}
    inf_low = -10^10;
    inf_up = 10^10;
    constrainttolerance = 10^-15;
    maxiter = 3000;
    maxfunevals = 3000;
    flag_display = false;
    fourier_normalizer_const()
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
        if strcmp(condi_Diri,'upbound') || strcmp(condi_Diri,'lowbound')  %'Dirichlet' only! %'Dirichlet' only!
            A2d = zeros(M,M);
            for ii = 1:M
                A2d(ii,:) = sin(2*pi*(M:N-1)*(ii-1)/N);
            end
            A_opt = (J_M.*Phi_M_c+J_M.^2.*Phi_M_s*bopi/(e-s) - J_M.^2.*Phi_M_smpn*bopi/(e-s));
            alt_Mm1 = ones(1,M);
            alt_Mm1(2:2:end)=-1;
            A_opt =A_opt.*alt_Mm1; 
            A_opt = (4/N)*(A_opt*A2d);
            if strcmp(condi_Diri,'upbound')
                A_opt = -A_opt;
                %b_opt = -((greeks(2)-greeks(1))/(e-s)-cond_b)/bopi;
                b_opt = -((greeks(2)-greeks(1))/(e-s)-up_val)/bopi;
            else
                %b_opt =((greeks(2)-greeks(1))/(e-s)-cond_b)/bopi;
                b_opt =((greeks(2)-greeks(1))/(e-s)-low_val)/bopi;
            end
        end
    end
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
    a_0_cons = 1/delta*(m_22*greeks(1) - m_12*greeks(2));
    a_1_cons = 1/delta*(-m_21*greeks(1) + m_11*greeks(2));
    zlast = [];
    myobjlast = [];
    glast = [];
    hlast = [];
    dmyobjlast = [];
    dglast = [];
    dhlast = [];
    
    if condi_v_low ||  condi_v_up

        b_opt = up_val*ones(M,1);
        A_opt = zeros(M,M);
        for kk =M+1:N
            for tt = M+1:N
                A_opt(kk-M,tt-M) = a_1_grad(tt-M)+ x_R(kk-M)*a_0_grad(tt-M) -4*b^2/(pi^2*N)*sum(J_M.^2.* sin(2*pi*(0:M-1)*(kk-1)/N).*sin(2*pi*(0:M-1)*(tt-1)/N)); 
            end
            b_opt(kk-M) = b_opt(kk-M) - a_1_cons - a_0_cons*x_R(kk-M);
        end      
        if condi_on_grid
            b_opt(1:m) = inf_up;
            b_opt(m+n+2:end) = inf_up;
        end
        if condi_v_low
            A_opt = -A_opt;
            b_opt = -low_val*ones(M,1);
            for kk =M+1:N
                b_opt(kk-M) = b_opt(kk-M) + a_1_cons + a_0_cons*x_R(kk-M);
            end
            if condi_on_grid
                b_opt(1:m) = inf_up;
                b_opt(m+n+2:end) = inf_up;
            end
        end
    end
    
    if condi_u_low ||  condi_u_up
        A_opt = zeros(M,M);
        b_opt = up_val*ones(M,1);
        for kk =M+1:N
            for tt = M+1:N
                A_opt(kk-M,tt-M) = a_0_grad(tt-M) -4*b/(pi*N)*sum(J_M.* cos(2*pi*(0:M-1)*(kk-1)/N).*sin( 2*pi*(0:M-1)*(tt-1)/N)) ; 
            end
            b_opt(kk-M) = b_opt(kk-M)-a_0_cons;
        end      
        if condi_on_grid
            b_opt(1:m) = inf_up;
            b_opt(m+n+2:end) = inf_up;
        end
        if condi_u_low
            A_opt = -A_opt;
            b_opt = -low_val*ones(M,1);
            for kk =M+1:N
                b_opt(kk-M) = b_opt(kk-M) + a_0_cons;
            end
            if condi_on_grid
                b_opt(1:m) = inf_up;
                b_opt(m+n+2:end) = inf_up;            
            end
        end
    end
    lb = [];
    ub = [];
    if condi_z_low
        lb = low_val*ones(M,1);
        if condi_on_grid
            lb(m+n+2:end) = inf_low;
            lb(1:m) = inf_low;
        end
    end
    if condi_z_up;
        ub = up_val*ones(M,1);  
        if condi_on_grid
            ub(m+n+2:end) = inf_up;
            ub(1:m) = inf_up;
        end
    end
    Aeq_opt = [];
    beq_opt = [];
    
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
   
   
    %counter = 1;
    %record = zeros(maxiter,1);
    
    function [myobj,g,h,dmyobj,dg,dh] = objcon(z_M) 
        z_N = [0,-fliplr(z_M(2:M)), z_M];  
        ifft_this = imag(ifft(z_N));
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
        F_M = fun(xvu_M, r_values); 
        myobj = sum((z_M-F_M).^2)/N;
        myobj = myobj';
        %handle u=y'
        DF_u_M = partial_u(xvu_M);
        Psi_u_M = (z_M-F_M).*DF_u_M;
        Psi_u_N = zeros(1,N);
        Psi_u_N(M+1:N) = Psi_u_M;
        I_u = sum (Psi_u_M);
        phi_u = 4*N*bopi*imag(ifft(J_N.*real(ifft(Psi_u_N))));
        phi_u = I_u*a_0_grad - phi_u(M+1:N);
        %handle v=y
        DF_v_M = partial_v(xvu_M);
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
    try
        [opt, myobjopt_inte,~, ~] = fmincon(@obj, z_init, A_opt,b_opt, Aeq_opt, beq_opt, lb, ub, @con, options);
    catch
        myobjopt_inte = -1;
    end
    if myobjopt_inte~= -1
        X_inte = [x_R; v_R; u_R; opt];
    else
        X_inte = [x_R; v_R; u_R; z_init];
    end
    %{
    temp = A_opt*opt' - b_opt;
    temp2 = A_opt*opt' + a_0_cons*x_R'  + a_1_cons*ones(M,1);
    z_N = [0,-fliplr(opt(2:M)), opt];
    temp3 = imag(ifft(J_N.^2.*imag(ifft(z_N))));
    temp4 = a_1*ones(1,M)+a_0*x_R - 2*b^2*N/(pi^2)*temp3(M+1:N);
    temp4 - v_R;
    A_check = zeros(M,M);
    for kk =M+1:N
        for tt = M+1:N
            A_check(kk-M,tt-M) = -4*b^2/(pi^2*N)*sum(J_M.^2.*sin(2*pi*(0:M-1)*(kk-1)/N).*sin(2*pi*(0:M-1)*(tt-1)/N)); 
        end
    end   
    temp5 = A_check*opt' + a_0*x_R' +a_1*ones(M,1);
    %}
end

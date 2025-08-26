%ode_order_1_by_z_engine.m
%{
implementation details with opt on z.  see doc
Implementation specification: optimization on $z$ ODE system of order 1
%}

function [u_N_out, fval] = ode_order_1_by_z_engine(fun, partial_u, get_H, get_H_der1,...
    y_s_vect, s, e, p, init_guess, ws, pps, flag_applying_constrain, weight_constrain, r_val, test_type, ps, qs, c3, M, co_R)
    global options
    fourier_normalizer_const()
    [M, dim] = size(init_guess);
    N = 2*M;
    n = 2^p;
    lambda = (e-s)/n;
    m = (M-n)/2;
    delta = lambda*m;
    o = s - delta;
    b = e-s + 2*delta;
    x_N = -b + o + lambda*(0:N-1);
    x_R = x_N(M+1:end);
    www = zeros(M,1);
    if flag_applying_constrain
        www(m+1:m+n+1) = 1;  %used for H term (constrain) 
    end
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
    opt_weight = 1;
    function [myobj_comb,g,h,dmyobj_comb,dg,dh] = objcon(z_comb)
        dmyobj_comb = zeros(dim*M,1);
        myobj_comb = 0;
        xu_comb_right = zeros(M, dim+1); % for the righ half !!!
        xu_comb_right(:,1) = x_R;
        for alpha = 1:dim
            z_M = z_comb((M*(alpha-1)+1):M*alpha);
            z_N = [z_M,0, -fliplr(z_M(2:M))];     
            a_0 = y_s_vect(alpha) + 2*bopi * dot(Phi_N, imag(ifft(z_N))); % Eq. (15)
            u_N = a_0 - 2*bopi*N*real(ifft(J_N.*(imag(ifft(z_N)))));  % Eq. (14)
            xu_comb_right(:,alpha+1) = u_N(M+1:end)';
        end
        
        F_M_comb_right = fun(xu_comb_right,r_val, test_type, ps, qs, M, co_R); 
        if flag_applying_constrain
            H_right = get_H(xu_comb_right(:,2:dim+1), test_type, c3); %new for constrain 
            H_left = [0;fliplr(H_right(2:M)')'];
            H_left = H_left.*www;
        end
        for alpha = 1:dim %index of u, i.e. alpha in doc
            if flag_applying_constrain
                DH_alpha_right = get_H_der1(xu_comb_right(:,2:dim+1), alpha, test_type, c3); 
                constrain_comp_alpha_right = weight_constrain*H_right.*DH_alpha_right; 
                constrain_comp_alpha_left = [0;fliplr(constrain_comp_alpha_right(2:M)')']; 
                constrain_comp_alpha_left = constrain_comp_alpha_left.*www;
            end
            Psi_M_alpha_left = zeros(M,1);
            for beta = 1:dim %index of F, i.e. beta in doc
                DF_M_right = partial_u(xu_comb_right, beta, alpha,test_type, co_R, ps, qs);
                DF_M_left = -[0;fliplr(DF_M_right(2:M)')'];
                z_M_left = z_comb((M*(beta-1)+1):M*beta);
                F_M_right = F_M_comb_right(:,beta);
                F_M_left = -[0;fliplr(F_M_right(2:M)')'];
                Psi_M_alpha_left = Psi_M_alpha_left + ws(beta)*(pps(beta)*z_M_left'-F_M_left).*DF_M_left;  %by definition of Phi_alpha
            end
            if flag_applying_constrain
                Psi_M_alpha_left = Psi_M_alpha_left -constrain_comp_alpha_left; %by definition of Phi_alpha
            end
            I_alpha =  sum(Psi_M_alpha_left); %by definition of I_alpha
            Psi_N_alpha = zeros(1,N);
            Psi_N_alpha(1:M) = Psi_M_alpha_left;
            F_M_right = F_M_comb_right(:,alpha);
            F_M_left = -[0;fliplr(F_M_right(2:M)')'];
            z_M_left = z_comb((M*(alpha-1)+1):M*alpha);
            myobj_comb = myobj_comb + sum((z_M_left-F_M_left').^2)/(dim*N);
            if flag_applying_constrain
                myobj_comb = myobj_comb  + weight_constrain*sum(H_left.^2)/(dim*N); 
            end            
            W = 4*bopi*I_alpha*ifft_Phi - 4*bopi*N*imag(ifft(J_N.*real(ifft(Psi_N_alpha)))); %Eq. (18)
            dmyobj_comb((M*(alpha-1)+1):M*alpha)= (ws(alpha)*pps(alpha)*(pps(alpha)*z_M_left-F_M_left')-W(1:M))/(dim*M); %Eq. (19-20)
        end
        myobj_comb = myobj_comb*opt_weight;
        dmyobj_comb = dmyobj_comb*opt_weight;

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
    z_init = zeros(1, dim*M);
    for alpha_ = 1:dim
        z_init((M*(alpha_-1)+1):M*alpha_) = init_guess(:,alpha_);
    end
    [zopt, fval] = fmincon(@obj, z_init, [],[], [], [], lb, ub, @con, options);
    u_N_out = zeros(N, dim);
    for a = 1:dim
        z = zopt((M*(a-1)+1):M*a);
        z = [z,0, -fliplr(z(2:M))];     
        a_0_this = y_s_vect(a) + 2*bopi * dot(Phi_N, imag(ifft(z)));
        u = a_0_this - 2*bopi*N*real(ifft(J_N.*(imag(ifft(z))))); 
        u_N_out(:, a) = u;
    end
end

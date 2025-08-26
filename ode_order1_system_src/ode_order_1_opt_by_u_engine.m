%ode_order_1_opt_by_u_engine.m
%{
Engine for TIBO on non-linear ODE system. 
See following model specification document:
"Implementation specification: optimization on $u$ ODE system of order $1$"
%}

function [u_tibo, fval] = ode_order_1_opt_by_u_engine(fun, partial_u, get_H, get_H_der1,...
    y_s_vect, s, e, p, init_guess, flag_applying_constrain, weight_constrain, r_val, test_type, ps, qs, c3, M, co_R)
    global options
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
    R=ones(M,1);
    %R = zeros(M,1);
    %R(m+1:m+n+1)=1;
    www = zeros(M,1);
    if flag_applying_constrain
        www(m+1:m+n+1) = 1;  %used for H term (constrain) 
    end
    C_M = -4*ones(M,1);
    C_M(1) = -2;
    J_M = 0:(M-1);
    J_N = zeros(1,N);
    J_N(1:M) = J_M;
    piob = pi/b;
    
    zlast = [];
    myobjlast = [];
    glast = [];
    hlast = [];
    dmyobjlast = [];
    dglast = [];
    dhlast = [];
    lb = [];
    ub = [];
    O_M = zeros(M,1);

    function [myobj_comb,g,h,dmyobj_comb,dg,dh] = objcon(u_comb)
        dmyobj_comb = zeros(dim*(M-1),1);
        myobj_comb = 0;
        xu_comb_right = zeros(M, dim+1); % for the righ half, covering [0,b]
        xu_comb_right(:,1) = x_R;
        z_combin_right = zeros(M, dim);
        for alpha = 1:dim
            u_R = u_comb(((M-1)*(alpha-1)+1):(M-1)*alpha);
            u_R = [u_R(1:m), y_s_vect(alpha),  u_R((m+1):end)];
            u_N = [0,fliplr(u_R(2:M)), u_R(1:M)];     
            z_N =-2*pi*N/b*imag(ifft(real(ifft(u_N)).*J_N));  %Eq (15)
            z_combin_right(1:M,alpha) = z_N((M+1):end);
            xu_comb_right(:,alpha+1) = u_R';
        end
        
        F_comb_right = fun(xu_comb_right,r_val,test_type, ps, qs, M, co_R);
        if flag_applying_constrain
            H_right = get_H(xu_comb_right(:,2:(dim+1)), test_type, c3);
        end
        for alpha = 1:dim %index of u, i.e. alpha in doc
            if flag_applying_constrain
                DH_alpha_right = get_H_der1(xu_comb_right(:,2:(dim+1)), alpha,test_type, c3); %new for constrain 
                constrain_comp_alpha_right = weight_constrain*H_right.*DH_alpha_right.*www; %new for constrain 
            end
            Psi_alpha_right = zeros(M,1);
            for beta = 1:dim %index of F, i.e. beta in doc
                DF_R = partial_u(xu_comb_right, beta, alpha,test_type, co_R, ps, qs);
                z_R = z_combin_right(:,beta);
                F_R = F_comb_right(:,beta);
                Psi_alpha_right = Psi_alpha_right + R.*(z_R-F_R).*DF_R;   %by Definition of Psi_{alpha}
            end
            if flag_applying_constrain
                Psi_alpha_right = Psi_alpha_right -constrain_comp_alpha_right; 
            end
            F_right = F_comb_right(:,alpha);
            z_right = z_combin_right(:,alpha);
            z_F = (z_right-F_right).*R;
            
            W = N*piob*real(ifft(imag(ifft([O_M; z_F])).*J_N'));
            W = W((M+1):N).*C_M;  %Eq (19)
            grad_part_alpha = (W - Psi_alpha_right)/2;  %Eq. (20)
            grad_part_alpha = [grad_part_alpha(1:m); grad_part_alpha((m+2):end)];
            dmyobj_comb(((M-1)*(alpha-1)+1):(M-1)*alpha)= grad_part_alpha; % Eq. (21)
            
            z_F = [z_F(1:m); z_F((m+2):end)];
            myobj_comb = myobj_comb + sum(z_F.^2);
            %sprintf('max error %0.6e alpha: %d\n', max(abs(z_F)), alpha)
            if flag_applying_constrain
                H_right_minus_s = [H_right(1:m); H_right((m+2):end)];
                myobj_comb = myobj_comb  + weight_constrain*sum(H_right_minus_s.^2);
            end            
        end
        %myobj_comb = myobj_comb';
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

    u_init = zeros(1, dim*(M-1));
    init_guess = [init_guess(1:m,:); init_guess(m+2:end,:)]; %remove rows associated to s    
    for alphat = 1:dim
        u_init(((M-1)*(alphat-1)+1):(M-1)*alphat) = init_guess(:,alphat);
    end
    [uopt, fval] = fmincon(@obj, u_init, [],[], [], [], lb, ub, @con, options);
    u_tibo = zeros(M, dim);
    for a = 1:dim
        u = uopt(((M-1)*(a-1)+1):(M-1)*a);
        u = [u(1:m), y_s_vect(a),  u((m+1):end)];
        u_tibo(:, a) = u;
    end
end

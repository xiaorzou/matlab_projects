%fourier_normalizer_app_ode_deg_1_system_v0.m
%{
this is the version for paper data
see model specification document
%}

function [u_tibo, fval] = ode_order_1_by_u_enhance_engine(fun, partial_u, get_H, get_H_der1,...
    y_s_vect, s, e, p, init_guess, ws, pps, flag_applying_constrain, weight_constrain, r_val, test_type, ps, qs, c3, M, co_R)
    %addpath('D:/matlab/fourier_normalizer/fourier_normalizer_engine')  
    %addpath('D:/matlab/cos_approx/cos_approx_engine')  
    global options
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
    x_R_b = [x_R, b];
    www = zeros(M+1,1);
    %pos_right = (N-(m+n+1)+2): (N-(m+1)+2);
    pos_left = m+1:m+n+1;
    if flag_applying_constrain
        www(pos_left) = 1;  %used for H term (constrain) 
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
    ifft_Phi = ifft_Phi(1:M+1);
    
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
    
    % u_comb is for the left part:  u_comb(1)..., u_comb(M) are assoicte to
    % the first half of u vector with N element), i.e. they are assoicated
    % to x_1=-b, ..., x_M=-lambda, note that x_{M+1}=0. 
    %function [myobj_comb,g,h,dmyobj_comb,dg,dh] = objcon(z_comb)
    function [myobj_comb,g,h,dmyobj_comb,dg,dh] = objcon(u_comb)
        dmyobj_comb = zeros(dim*M,1);
        myobj_comb = 0;
        xu_comb_right = zeros(M+1, dim+1); % for the righ half !!!
        xu_comb_right(:,1) = x_R_b;
        for alpha = 1:dim
            %z_M = z_comb((M*(alpha-1)+1):M*alpha);
            %z_N = [z_M,0, -fliplr(z_M(2:M))];     
            %a_0 = y_s_vect(alpha) + 2*bopi * dot(Phi_N, imag(ifft(z_N)));
            %u_N = a_0 - 2*bopi*N*real(ifft(J_N.*(imag(ifft(z_N))))); % we need u_N, rather than just u_M, u(0) is not known!
            u_M = u_comb((M*(alpha-1)+1):M*alpha);
            u_M_b = [u_M(1:m+n), y_s_vect(alpha), u_M(m+n+1:end)];
            u_M_b = fliplr(u_M_b);
            xu_comb_right(:,alpha+1) = u_M_b;
        end
        u_left = (fliplr(xu_comb_right(:, 2:end)'))';
        F_M_comb_right = zeros(M+1,dim);
        F_M_comb_right(1:M,:) = fun(xu_comb_right(1:M,:),r_val, test_type, ps, qs, M, co_R);  %the values on the right half!!!
        if flag_applying_constrain
            H_right = get_H(xu_comb_right(:,2:dim+1), test_type, c3); %new for constrain 
            H_left = fliplr(H_right')';
        end
        u_tilda_comb = zeros(M+1, dim);
        for beta = 1:dim
            F_beta = -fliplr(F_M_comb_right(:,beta)');
            a_0 = y_s_vect(beta)+ (4*b/pi)*dot(ifft_Phi, F_beta); 
            temp = (real(ifft(imag(ifft([F_beta, zeros(1,M-1)])).*J_N)))'; 
            u_tilda_comb(:,beta) = a_0*ones(M+1,1) - (4*b*N/pi)*  temp(1:M+1);
        end
        u_delta_comb = u_tilda_comb - u_left;
        s_vect = sum(u_delta_comb);
        
        for alpha = 1:dim %index of u, i.e. alpha in doc
            if flag_applying_constrain
                DH_alpha_right = get_H_der1(xu_comb_right(:,2:dim+1), alpha, test_type, c3); %new for constrain 
                DH_alpha_left = fliplr(DH_alpha_right')';
                constrain_comp_alpha_left = H_left.*DH_alpha_left.*www;
            end
            %Psi_M_alpha_left = zeros(M,1);
            eta_alpha = zeros(M+1,1);
            for beta = 1:dim %index of F, i.e. beta in doc
                DF_beta_alpha = partial_u(xu_comb_right(1:M,:), beta, alpha,test_type, co_R, ps, qs);
                DF_beta_alpha = -[0;fliplr(DF_beta_alpha')'];
                part_double_ifft = imag(ifft(real(ifft([u_delta_comb(:, beta); zeros(M-1,1)])).*J_N'));
                part_double_ifft = part_double_ifft(1:M+1);
                a_beta_0_grad =  (4*b/pi)* ifft_Phi'.*DF_beta_alpha;
                eta_alpha = eta_alpha + s_vect(beta)*a_beta_0_grad - (4*b*N/pi)*DF_beta_alpha.*part_double_ifft;   
            end
            gradient_alpha = eta_alpha;
            if flag_applying_constrain
                gradient_alpha = gradient_alpha + constrain_comp_alpha_left; %new for constrain 
            end
            
            myobj_comb = myobj_comb + sum((u_delta_comb(:,alpha)).^2);
            %sprintf('max error %0.6e alpha: %d\n', max(abs(z_M_left-F_M_left')), alpha)
            if flag_applying_constrain
                myobj_comb = myobj_comb  + weight_constrain*sum(H_left.^2); %new for constrain 
            end
           
            gradient_alpha = gradient_alpha -u_delta_comb(:,alpha);
            gradient_alpha = [gradient_alpha(1:m+n); gradient_alpha(m+n+2:end)];
            
            dmyobj_comb((M*(alpha-1)+1):M*alpha)= gradient_alpha;
        end
        myobj_comb = myobj_comb*opt_weight;
        dmyobj_comb = dmyobj_comb*opt_weight;
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
    %z_init = zeros(1, dim*M);
    u_init = zeros(1, dim*M);
    for alpha_ = 1:dim
        u_init((M*(alpha_-1)+1):M*alpha_) = init_guess(:,alpha_);
    end
    [uopt, fval] = fmincon(@obj, u_init, [],[], [], [], lb, ub, @con, options);
    u_tibo = zeros(M+1, dim);
    for a = 1:dim
        u = uopt((M*(a-1)+1):M*a);
        u_tibo(:,a) = [u(1:m+n), y_s_vect(a),  u((m+n+1):end)];
    end
end

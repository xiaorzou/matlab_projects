
%{
issue to be clarifed

for mix,
myobjopt_inte can only to 10^-9 with default alpha and beta

if beta = beta_detaul +1,  goest to 10^-17,  but diff is still 10^-3?

%}

function V = fourier_normalizer_app_ode_order2_en_linear(X_inte, a_0, a_1, b, G, Q, R, alpha, beta, m, n, cond_type,  A1,A2)
    Q = [Q,0];
    R = [R,0];
    z_M = X_inte(4,:);
    M = length(z_M);
    N = 2*M;
    z_N = [0,-fliplr(z_M(2:M)), z_M];  
    J_M = zeros(1,M);
    J_M(2:M) = 1./(1:M-1);
    A = ones(1,N);
    A(2:2:end) = -1;
    b_M = 2*imag(ifft(z_N));
    b_M = b_M.*A;
    b_M = b_M(1:M);
    b_0 = b_M(1);
    b_M = b_M(2:end);
    x = X_inte(1,:);    
    x_R_plusb = zeros(1, M+1);
    x_R_plusb(1:M) = x;
    x_R_plusb(M+1) = b;
    x_0 = x(1);
    x = x(2:end);
    v = X_inte(2,:);
    v_Mplus1 = zeros(1,M+1);
    u_Mplus1 = zeros(1,M+1);
    v_Mplus1(1:M) = v;
    
    v_0 = v(1);
    v_M = a_1 + a_0*b;
    v_Mplus1(M+1) = v_M;
    v = v(2:end);
    u = X_inte(3,:);
    u_Mplus1(1:M) = u;
    u_0 = u(1);
    u = u(2:end);
    z = X_inte(4,:);
    z_0 = z(1);
    z = z(2:end);
    bopi = b/pi;
    A = zeros(M-1,M-1);
    C = zeros(M-1,M-1);
    K = 1:(M-1);
    I = ones(M-1,1);
    alt = I';
    alt(1:2:end)=-1;
    alt = -alt;
    
    %v_M = a_1+b*a_0;
    u_M = a_0 - bopi*(-1)*sum(alt.*b_M.*(1./K));
    u_Mplus1(end) = u_M;
    for i = 1:M-1
        A(i,:) = sin(2*pi*K*i/N);
        C(i,:) = cos(2*pi*K*i/N);
    end
    
    lambda = 2*b/N;
    b_check = 2*a_1/(bopi^2*M)*diag(K.^2)*A*I + 2*lambda*a_0/(M*bopi^2)*diag(K.^2)*A*K' - 2/(M*bopi^2)*diag(K.^2)*A*v';
    b_check_result = b_M'-b_check;
    u_check = v_0*(2/(M^2*bopi)*C*diag(K)*A*K' - I/b-2/(M*bopi)*C*diag(K)*A*I) + v_M*(I/b - 2/(M^2*bopi) *C *diag(K)*A*K') + 2/(bopi*M)*C*diag(K)*A*v';
    u_check_result = u'-u_check;
    

    u_0_check_1 =v_0*(-1/b + 1/(M*bopi)*sum(alt.*K.*cot(pi*K/N)) - 1/bopi*sum(alt.*cot(pi*K/N)));
    u_0_check_2 =  v_M *(1/b - 1/(M*bopi)*sum(alt.*K.*cot(pi*K/N)));
    u_0_check_3 =  1/(bopi)*sum(alt.*v.*cot(K*pi/N));
    

     u_check_0_2_1 = v_0*(2/(M^2*bopi)*I'*diag(K)*A*K' - 1/b-2/(M*bopi)*I'*diag(K)*A*I); 
     u_check_0_2_2 = v_M*(1/b - 2/(M^2*bopi) *I' *diag(K)*A*K');
     u_check_0_2_3 =  2/(bopi*M)*I'*diag(K)*A*v';
    r1 = u_0_check_1 - u_check_0_2_1;
    r2 = u_0_check_2 - u_check_0_2_2;
    r2 = u_0_check_3 - u_check_0_2_3;
    u_0_lim = (v_0-v_M)/(1-2/b);
    
    %v_M = a_1+b*a_0;
    %u_M = a_0 - bopi*(-1)*sum(alt.*b_M.*(1./K));
    
    u_check_M_2_1 = v_0*(2/(M^2*bopi)*(-alt)*diag(K)*A*K' - 1/b-2/(M*bopi)*(-alt)*diag(K)*A*I); 
    u_check_M_2_2 = v_M*(1/b - 2/(M^2*bopi) *(-alt) *diag(K)*A*K');
    u_check_M_2_3 =  2/(bopi*M)*(-alt)*diag(K)*A*v';
    u_M_r = u_check_M_2_1 + u_check_M_2_2 + u_check_M_2_3;
    
    check_u_M_1 = -v_0*(1/b+1/(M*bopi)*sum(K.*cot(K*pi/N))) +v_0*1/(bopi)*sum((alt).*tan(K*pi/N));
    check_u_M_2 = v_M*(1/b+1/(M*bopi)*sum(K.*cot(K*pi/N)));
    check_u_M_3 = -1/(bopi)*sum((alt).*v.*tan(K*pi/N));
    check_u_M_r = check_u_M_1 + check_u_M_2 + check_u_M_3;
    
    u_check_details = zeros(M-1,1);
    for i = 1:M-1
        cot_i = cot((K+i)*pi/N) + cot((K-i)*pi/N);
        cot_i(i) = cot((i+i)*pi/N);
         temp1 = v_0 *(0.5/bopi*sum((-1)^(i+1)*alt.*cot_i) - 0.5/(bopi*M)*sum((-1)^(i+1)*alt.*cot_i.*K) -1/b);
         temp2 = v_M*(0.5/(bopi*M)*sum((-1)^(i+1)*alt.*cot_i.*K)  +1/b );
         temp3 = -0.5/bopi*sum((-1)^(i+1)*alt.*cot_i.*v);
         u_check_details(i,1) = temp1 +temp2 +temp3;
    end
    disp(max(abs(u_check - u_check_details)));
    % U = Alpha*V,  
    Alpha = zeros(1+M, 1+M);
    Alpha(1,1)= 1/(M*bopi)*sum(alt.*K.*cot(K*pi/N))-1/b-1/bopi*sum(alt.*cot(K*pi/N));
    
    Alpha(1,M+1) = -1/(M*bopi)*sum(alt.*K.*cot(K*pi/N))+1/b;
    Alpha(1,2:M) = 1/bopi*alt.*cot(K*pi/N);
    Alpha(M+1,1) = -1/(M*bopi)*sum(K.*cot(K*pi/N))-1/b+1/bopi*sum(alt.*tan(K*pi/N));
    Alpha(M+1,M+1) = 1/b+1/(M*bopi)*sum(K.*cot(K*pi/N));
    Alpha(M+1,2:M) = -1/bopi*alt.*tan(K*pi/N);
    for i = 1:M-1
        cot_i = cot((K+i)*pi/N) + cot((K-i)*pi/N);
        cot_i(i) = cot((i+i)*pi/N);
        Alpha(i+1,1) = 0.5/bopi*sum((-1)^(i+1)*alt.*cot_i) - 0.5/(bopi*M)*sum((-1)^(i+1)*alt.*cot_i.*K) -1/b;
        Alpha(i+1,M+1) = 0.5/(bopi*M)*sum((-1)^(i+1)*alt.*cot_i.*K)  +1/b ; 
        Alpha(i+1,2:M) = -0.5/bopi*(-1)^(i+1)*alt.*cot_i;
    end
    disp(max(abs(u_Mplus1' - Alpha*v_Mplus1')));
    
    O = (2/M)^0.5*A;
    Theta = O*diag(1./(K.^2))*O;

    
    g = zeros(1,M+1);
    g(1) = -alpha;
    g(2:M) = Theta*G(2:M)';
    g(M+1) = -beta;
    g = -g;
    
    
    Phi = zeros(M+1,M+1);
    %{ 
    Phi(1,:) = 0;
    Phi(1,m+1) = 1;
    Phi(1+M,:) = 0;
    Phi(1+M,m+n+1) = 1;
    %}

    if strcmp(cond_type, 'Dirichlet')
        Phi(1,:) = 0;
        Phi(1,m+1) = 1;
        Phi(1+M,:) = 0;
        Phi(1+M,m+n+1) = 1;
    else
        Phi(1,:) = A1(2)*Alpha(m+1,:) + A1(4)*Alpha(m+n+1,:);
        Phi(1,m+1) = Phi(1,m+1) +  A1(1);
        Phi(1,m+n+1) = Phi(1,m+n+1) +  A1(3);
        Phi(M+1,:) = A2(2)*Alpha(m+1,:) + A2(4)*Alpha(m+n+1,:);
        Phi(M+1,m+1) = Phi(M+1,m+1) +  A2(1);
        Phi(M+1,m+n+1) = Phi(M+1,m+n+1) +  A2(3);
        % first row assolucat to condition 1
    end
    
    
    for i = 1:(M-1)
        Phi(i+1,1) = -(M-i)/(M*bopi^2) + (Theta(i,:).*R(2:M))*Alpha(2:M,1);
        Phi(i+1,1+M) = -i/(M*bopi^2) + (Theta(i,:).*R(2:M))*Alpha(2:M,1+M);
        Phi(i+1,2:M) = Theta(i,:).*Q(2:M) +  (Theta(i,:).*R(2:M))*Alpha(2:M,2:M);
        Phi(i+1,i+1) = Phi(i+1,i+1) + 1/bopi^2;
    end
    %V = inv(Phi)*g';
    V = Phi\g';
    %disp(max(abs(V-v_Mplus1')));
    %pos = find(x_R_plusb>=1 & x_R_plusb<=3);
    %V_res = V(pos);
end
%{
det(Phi) = 2.551483076220252e-38
%}
    
    

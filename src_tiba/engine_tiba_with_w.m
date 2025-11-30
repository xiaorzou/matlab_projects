
%{
the engine for the following paper
"An Application of the Trigonometric Interpolation on Second Order Linear ODE"
%}

function V = engine_tiba_with_w(b, W, R, Q, P, alpha, beta, m, n, cond_type,  A1,A2, v_s)
    %Q = [Q,0];
    %P = [P,0];
    M = length(R);
    N = 2*M;
    bopi = b/pi;   
    S = zeros(M-1,M-1);
    C = zeros(M-1,M-1);
    K = 1:(M-1);
    I = ones(M-1,1);
    alt = I';
    alt(1:2:end)=-1;
    alt = -alt;
    for i = 1:M-1
        S(i,:) = sin(2*pi*K*i/N);
        C(i,:) = cos(2*pi*K*i/N);
    end
    A = zeros(1+M, 1+M);
    A(1,1)= 1/(M*bopi)*sum(alt.*K.*cot(K*pi/N))-1/b-1/bopi*sum(alt.*cot(K*pi/N));
    
    A(1,M+1) = -1/(M*bopi)*sum(alt.*K.*cot(K*pi/N))+1/b;
    A(1,2:M) = 1/bopi*alt.*cot(K*pi/N);
    A(M+1,1) = -1/(M*bopi)*sum(K.*cot(K*pi/N))-1/b+1/bopi*sum(alt.*tan(K*pi/N));
    A(M+1,M+1) = 1/b+1/(M*bopi)*sum(K.*cot(K*pi/N));
    A(M+1,2:M) = -1/bopi*alt.*tan(K*pi/N);
    for i = 1:M-1
        cot_i = cot((K+i)*pi/N) + cot((K-i)*pi/N);
        cot_i(i) = cot((i+i)*pi/N);
        A(i+1,1) = 0.5/bopi*sum((-1)^(i+1)*alt.*cot_i) - 0.5/(bopi*M)*sum((-1)^(i+1)*alt.*cot_i.*K) -1/b;
        A(i+1,M+1) = 0.5/(bopi*M)*sum((-1)^(i+1)*alt.*cot_i.*K)  +1/b ; 
        A(i+1,2:M) = -0.5/bopi*(-1)^(i+1)*alt.*cot_i;
    end
    
    O = (2/M)^0.5*S;
    %Theta = O*diag(1./(K.^2))*O;
    Theta_inv = O*diag(K.^2)*O;
    
    Psi = zeros(1,M+1);
    Psi(1) = -alpha;
    %Psi(2:M) = Theta*R(2:M)'; %old
    Psi(2:M) = R(2:M)'; %new
    Psi(M+1) = -beta;
    Psi = -Psi;
   
    Phi = zeros(M+1,M+1);
    diag_W_hat = diag(W(2:end));
    diag_P_hat = diag(P(2:end));
    diag_Q_hat = diag(Q(2:end));

    A_hat = A(1:end-1,:);
    A_hat = A_hat(:,1:end-1);
    A_hat = A_hat(2:end,:);
    A_hat = A_hat(:,2:end);
    A_0 = A(2:end-1,1);
    A_M = A(2:end-1,end);
    C_hat =  (1/bopi^2)*diag_W_hat*Theta_inv + diag_Q_hat + diag_P_hat*A_hat;
    C_0 = diag_P_hat*A_0 - 1/(bopi^2*M)*diag_W_hat*Theta_inv*(M*I-K');
    C_M = diag_P_hat*A_M - 1/(bopi^2*M)*diag_W_hat*Theta_inv*K';
    
    %geneneral case should be consistent to Dirchilet
    if strcmp(cond_type, 'Dirichlet')
        Phi(1,:) = 0;
        Phi(1,m+1) = 1;
        Phi(1+M,:) = 0;
        Phi(1+M,m+n+1) = 1;
    else
        Phi(1,:) = A1(2)*A(m+1,:) + A1(4)*A(m+n+1,:);
        Phi(1,m+1) = Phi(1,m+1) +  A1(1);
        Phi(1,m+n+1) = Phi(1,m+n+1) +  A1(3);
        Phi(M+1,:) = A2(2)*A(m+1,:) + A2(4)*A(m+n+1,:);
        Phi(M+1,m+1) = Phi(M+1,m+1) +  A2(1);
        Phi(M+1,m+n+1) = Phi(M+1,m+n+1) +  A2(3);
    end
    Phi(2:end-1,2:end-1)=C_hat;
    Phi(2:end-1,1) = C_0;
    Phi(2:end-1,end) = C_M;
    
    
    zero_rows = [];
    for row = 1:(M+1)
        scaling = max(abs(Phi(row,:)));        
        if scaling ==0
            zero_rows = [zero_rows, row];
        else       
            Phi(row,:) = Phi(row,:)/scaling;
            Psi(row) = Psi(row)/scaling;
        end
    end
    
    %try to handle the situgation where Phi is singular when q is larger and co_R(2)=0
    
    zero_rows = fliplr(zero_rows);
    Phi_adj = Phi;
    for i = zero_rows
        if i==1
            Phi_adj(:,2) = Phi_adj(:,2) + Phi_adj(:,1); %assume that v_1=v_{2}
            Phi_adj = Phi_adj(2:end,:);
            Phi_adj = Phi_adj(:,2:end);
            Psi = Psi(2:end);
        elseif i == M+1
            Phi_adj(:,M) = Phi_adj(:,M) + Phi_adj(:,1+M); %assume that v_1=v_{2}
            Phi_adj = Phi_adj(1:end-1,:);
            Phi_adj = Phi_adj(:,1:end-1);
            Psi = Psi(1:end-1);
        elseif i<M/2
            Phi_adj(:,i+1) = Phi_adj(:,i+1) + Phi_adj(:,i); %assume that v_1=v_{2}
            Phi_adj = [Phi_adj(1:i-1,:); Phi_adj(i+1:end,:)];
            Phi_adj = [Phi_adj(:,1:i-1), Phi_adj(:,i+1:end)];
            Psi = [Psi(1:i-1), Psi(i+1:end)];
        else
            Phi_adj(:,i-1) = Phi_adj(:,i-1) + Phi_adj(:,i); %assume that v_1=v_{2}
            Phi_adj = [Phi_adj(1:i-1,:); Phi_adj(i+1:end,:)];
            Phi_adj = [Phi_adj(:,1:i-1), Phi_adj(:,i+1:end)];
            Psi = [Psi(1:i-1), Psi(i+1:end)];            
        end
    end
   
    Psi = Psi';
    %{  
    %does not solve sigularity
    if strcmp(cond_type, 'Neumann') && W(m+1)==0
        Psi = Psi - Phi_adj(:,m+1)*v_s;
        Psi = Psi(2:end,:);
        Phi_adj = Phi_adj(2:end,:);
        Phi_adj = [Phi_adj(:,1:m), Phi_adj(:,m+2:end)];
    end
    %}
    
    V = Phi_adj\Psi;
    %{
    %does not solve sigularity
    if strcmp(cond_type, 'Neumann') && W(m+1)==0
        V = [V(1:m); v_s; V(m+1:end)];
    end
    %}
   
    zero_rows = fliplr(zero_rows);
    for i = zero_rows
        if i==1
            V = [V(1); V];
        elseif i == M+1
            V = [V; V(end)];
        elseif i<M/2
            V = [V(1:i-1); V(i) ; V(i:end)];
        else
            V = [V(1:i-1); V(i-1) ; V(i:end)];
        end
    end
   
end

    
    

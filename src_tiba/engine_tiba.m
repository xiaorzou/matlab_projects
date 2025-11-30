
%{
the implementation of the algorithm described in the following paper
"Trigonometric Interpolation Based Based Approach for Second Order ODE with
Mixed Boundary Conditions"
%}

function V = engine_tiba(b, R, Q, P, alpha, beta, m, n, cond_type,  A1,A2)
    Q = [Q,0];
    P = [P,0];
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
    Theta = O*diag(1./(K.^2))*O;

    
    g = zeros(1,M+1);
    g(1) = -alpha;
    g(2:M) = Theta*R(2:M)';
    g(M+1) = -beta;
    g = -g;
   
    Phi = zeros(M+1,M+1);

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
    
    
    for i = 1:(M-1)
        Phi(i+1,1) = -(M-i)/(M*bopi^2) + (Theta(i,:).*P(2:M))*A(2:M,1);
        Phi(i+1,1+M) = -i/(M*bopi^2) + (Theta(i,:).*P(2:M))*A(2:M,1+M);
        Phi(i+1,2:M) = Theta(i,:).*Q(2:M) +  (Theta(i,:).*P(2:M))*A(2:M,2:M);
        Phi(i+1,i+1) = Phi(i+1,i+1) + 1/bopi^2;
    end
    V = Phi\g';
end

    

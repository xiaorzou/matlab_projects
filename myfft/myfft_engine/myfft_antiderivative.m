% q(x) = \sum_{0<=j<=N-1}\tilde(A)_j cos(pi*j*x/(right-left))
% A_j = (2/l)*real(\phi(j*pi)*exp(-j*pi*i*left/(right-left))),
% 0<=j<=N-1,
% \tilde(A)_j = A_j,  j>0,  \tilde(A)_0 = A_0/2
%

function FFTAD = myfft_antiderivative(density_coef, M_A, Left, Right)
%input:
%   @density_coef? same as \tilde(A) in above forulation
%   @M_A: max degree of antiderivatives
%   @Left/Right: the effective range of density q(x)
%output:
%   FFTAD: antiderivatives values upto degree M_D, at point left, left+lambda,
%   .... or D(1), ...D(N/2), notice 
%          a) D(N/4+1)=0 if left=-right
%          b) FFTAD(1,:)= cdf of q(x) 

    bma=Right-Left;
    bmapi=bma/pi;
    N=length(density_coef);
    Aones=ones(1,N);
    Aones(2:2:end)=-ones(1,N/2);
    K=0:1:N-1;
    Y=0:1:(N-1);
    Y=2*bma*Y/N+2*Left-Right; %same as Y*Lambda+2*Left-Right
    MA = ceil(M_A/2);
    FFTAD = zeros(MA*2,N);
    Im = zeros (MA, N);
    Re = zeros (MA, N);
    F1_new=density_coef(2:N)./K(2:N); %CDF
    F1_new=[0,F1_new];
    F1_new =bmapi*Aones.*F1_new;
    for i=1:MA
        Im(i,:) = N*imag(ifft(F1_new));
        F1_new(2:N) = bmapi*F1_new(2:N)./K(2:N);
        Re(i,:) = N*real(ifft(F1_new));
        F1_new(2:N) = bmapi*F1_new(2:N)./K(2:N);
    end
    factor = 1;
    for i=1:MA
        factor = factor *(2*i-1);
        FFTAD(2*i-1,:) = (-1)^(i-1)*Im(i,:) + (density_coef(1)/factor)* (Y-Left).^(i*2-1);
        factor_2 = 1;
        for m=1:(i-1)
            factor_2 = factor_2*(2*m-1); 
            FFTAD(2*i-1,:) = FFTAD(2*i-1,:) + (-1)^(i-m-1)/factor_2*(Y-Left).^(m*2-1)*Re(i-m, N/2+1);
            factor_2 = factor_2 * 2*m ;
        end   
        factor = factor * 2*i;
        FFTAD(2*i,:) = (-1)^(i)*Re(i,:) + (density_coef(1)/factor)* (Y-Left).^(i*2);
        factor_2 = 1;
        for m=1:i
            if m ~= 1
                factor_2 = factor_2*(2*m-2);
            end
            FFTAD(2*i,:) = FFTAD(2*i,:) + (-1)^(i-m)/factor_2*(Y-Left).^(m*2-2)*Re(i-m+1, N/2+1);
            factor_2 = factor_2 *(2*m-1);
        end       
    end
    FFTAD=FFTAD(:,N/2+1:N);
    FFTAD = FFTAD(1:M_A,:);
    FFTAD=max(FFTAD,0);
end


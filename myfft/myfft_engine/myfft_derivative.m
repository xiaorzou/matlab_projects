% q(x) = \sum_{0<=j<=N-1}\tilde(A)_j cos(pi*j*x/(right-left))
% A_j = (2/l)*real(\phi(j*pi)*exp(-j*pi*i*left/(right-left))),
% 0<=j<=N-1,
% \tilde(A)_j = A_j,  j>0,  \tilde(A)_0 = A_0/2

function FFTF = myfft_derivative(density_coef, M_D, Left, Right)
%input:
%   @density_coef? same as \tilde(A) in above forulation
%   @M_D: M_D-1 is the max degree of derivatives
%   @Left/Right: the effective range of density q(x)
%output:
%   FFTF: derivatives values upto degree M_D, at point left, left+lambda,
%   .... or D(1), ...D(N/2), notice 
%          a) D(N/4+1)=0 if left=-right
%          b) FFTF(1,:)=q(x) 
    N=length(density_coef);
    Aones=ones(1,N);
    Aones(2:2:end)=-ones(1,N/2);
    K=0:1:N-1;
    if mod(M_D,2)==0
        M = M_D+1;
    else
        M = M_D;
    end
    FFTF=zeros(M,N); % Derivative of q(x) up to order M-1, M has to be odd!
    F1=Aones.*density_coef;
    FFTF(1,:)=N*real(ifft(F1));  % order 0 derivative, q value
    for i=1:((M-1)/2)
        F1=F1.*K*pi/(Right-Left);  % 
        FFTF(2*i,:)=(-1)^i*N*imag(ifft(F1)); % order 2*i-1 derivative
        F1=F1.*K*pi/(Right-Left);  %  
        FFTF(2*i+1,:)=(-1)^i*N*real(ifft(F1)); % order 2*i derivative
    end
    FFTF=FFTF(:,N/2+1:N);
    FFTF=FFTF(1:M_D,:);
end


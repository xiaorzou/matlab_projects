% q(x) = \sum_{0<=j<=N-1}\tilde(A)_j cos(pi*j*x/(right-left))
% A_j = (2/l)*real(\phi(j*pi)*exp(-j*pi*i*left/(right-left))),
% 0<=j<=N-1,
% \tilde(A)_j = A_j,  j>0,  \tilde(A)_0 = A_0/2

function coefdensity=myfft_density_coef(model, inputVar, leftend,rightend, N)
%input
%   @model: 'Gauss', 'Merton',...
%   @inputVar: mapping providing model parameters 
%   @leftend/rightend? define effective range of density q(x)
%   @N: number of terms in cos expansion
%output
%   @coefdensity? same as \tilde(A) in above forulation


    x=0:1:(N-1);
    coefdensity=myfft_option_char(model, inputVar,  x*pi/(rightend-leftend)).*exp(-x*leftend*pi/(rightend-leftend)*complex(0,1));
    coefdensity=2*real (coefdensity)/(rightend-leftend);
    coefdensity(1)=coefdensity(1)/2;
end
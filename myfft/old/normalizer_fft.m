function val = normalizer_fft(alpha, beta, a, fftpower)
%{
test_case = 2
fftpower = 10;
a = 1;
alpha=8*pi;
beta = 16*pi;
N = 2^(fftpower);
L = 2*beta/N;
X = (0:1:(N-1))*L-beta;
if test_case == 1
    U = normalizer(X,alpha,beta, a);
elseif test_case == 2
    U = sin(X)+10*sin(10*X);
elseif test_case == 3
    U = X
elseif test_case == 4
    U = X.^3
end
Aones=ones(1,N);
Aones(2:2:end)=-ones(1,N/2);
B = Aones.* (2*imag(ifft(U)));
M=N/2; %ignore second half for possible bug, which can be captured by using above U =sin(X) + 10(10*X)
%M=N;
B = B(1:M);
Aones2=ones(1,M);
Aones2(2:2:end)=-ones(1,M/2);
L2 = 2*beta/M;
X2 = (0:1:(M-1))*L2-beta;
if test_case == 1
    U2 = normalizer(X2,alpha,beta, a);
elseif test_case == 2
    U2 = sin(X2)+10*sin(10*X2);
elseif test_case == 3
    U2 = X2
elseif test_case == 4
    U2 = X2.^3
end
B_alternative = B.*Aones2;
U_to_be_2 = M*imag(ifft(B_alternative));
diff = U2-U_to_be_2;
%}
N = 2^(fftpower);
L = 2*beta/N;
X = (0:1:(N-1))*L-beta;
U = normalizer(X,alpha,beta, a);
Aones=ones(1,N);
Aones(2:2:end)=-ones(1,N/2);
B = Aones.* (2*imag(ifft(U)));
%M=N/2; %ignore second half for possible bug, which can be captured by using above U =sin(X) + 10(10*X)
M = N;
B = B(1:M);
Aones2=ones(1,M);
Aones2(2:2:end)=-ones(1,M/2);
K = (beta/pi)./(1:1:M-1);
A0 = dot(K.*Aones2(2:end),B(2:end));
A = -K.*B(2:end).*Aones2(2:end);
A = [A0, A];
val = M*real(ifft(A));
val = val/max(val);
L2 = 2*beta/M;
X2 = (0:1:(M-1))*L2-beta;
val = [X2;val];
%plot(X2,val)
end

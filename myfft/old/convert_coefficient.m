function val = convert_coefficient(opt, fftpower,alpha,beta,a)

fftpower = 15;
a = 1;
alpha=15*pi;
beta = 16*pi;
opt = 1;

output_normalizer_fft = normalizer_fft(alpha, beta, a, fftpower);
X = output_normalizer_fft(1,:);
normalizer = output_normalizer_fft(2,:);
N = length(X);
f_X = test_function(X,opt);
U = normalizer.*f_X;
Aones=ones(1,N);
Aones(2:2:end)=-ones(1,N/2);
%if odd function!
B_N = Aones.* (2*imag(ifft(U)));
B = Aones.* (2*imag(ifft(f_X)));
M = N/2;
L2 = 2*beta/M;
X2 = (0:1:(M-1))*L2-beta;
X_small = X2(left_pos:right_pos);
left_pos = floor((beta-alpha)/L2);
right_pos = floor((beta+alpha)/L2);

B_N = B_N(1:M);
B = B(1:M);
Aones2=ones(1,M);
Aones2(2:2:end)=-ones(1,M/2);


B_N_alternative = B_N.*Aones2;
f_N_X_fft = M*imag(ifft(B_N_alternative));
f_X_2 = test_function(X2,opt);
f_X_2_small = f_X_2(left_pos:right_pos);
diff_N = f_X_2 -f_N_X_fft;
f_N_X_fft_small = f_N_X_fft(left_pos:right_pos);
diff_N_small = f_N_X_fft_small - f_X_2_small;

B_alternative = B.*Aones2;
f_X_fft = M*imag(ifft(B_alternative));
diff = f_X_2 -f_X_fft;
f_X_fft_small = f_X_fft(left_pos:right_pos);
diff_small = f_X_fft_small - f_X_2_small;


plot(X_small, f_N_X_fft_small, X_small, f_X_fft_small)












end


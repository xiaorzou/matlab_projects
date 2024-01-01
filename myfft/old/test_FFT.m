function output_args = test_FFT()
    F_coff = [1,1,0,0,0,0,1,1];
    fftpower = 4;
    N=2^fftpower;
    Left = 0;
    Right = pi; 
    M_A = 3;
    M_D = 3;
    F = [F_coff, zeros(1,N-length(F_coff))];
    Lambda=(Right-Left)/(N/2); %new
    D =(0:1:(N/2-1))*Lambda+Left;
    [D_fft, AD_fft] = FFT(F,N, M_A, M_D, Left, Right);
    output_args = [D_fft, AD_fft];
end
    



function output = test_fun(F, D, F_coff, order)
N = length(D);
deg = leng(F_coff);
output = zero(1,N);
    for i=1:N
        for j = 1:deg
            if order == 0
                output(i) = output(i) + F_coff(j)*cos((j-1)*D(i));
            elseif order == 1 %derivative order 1
                if j ~= 1
                    output(i) = output(i) - (j-1)*F_coff(j)*sin((j-1)*D(i));
                else
                    temp = 0;
                end
             elseif order == 2  %dervative order 2
                 if j ~= 1
                     output(i) = output(i) - (j-1)^2 *F_coff(j)*cos((j-1)*D(i));
                 end
             elseif order == 3  %dervative order 3
                 output(i) = output(i) + (j-1)^3 *F_coff(j)*sin((j-1)*D(i));
             elseif order == 4  %dervative order 4
                 output(i) = output(i) + (j-1)^4 *F_coff(j)*cos((j-1)*D(i));
             elseif order == -1 %derivative order 1
                 if j == 1
                     output(i) = output(i) + F_coff(j)* D(i)
                 else
                     output(i) = output(i) + (F_coff(j)/D(i))*sin((j-1)*D(i))
                 end
             elseif order == -2 %derivative order 1
                 if j == 1
                     output(i) = output(i) + F_coff(j)* D(i)^2/2
                 else
                     output(i) = output(i) - (F_coff(j)/D(i)^2)*cos((j-1)*D(i)) + (F_coff(j)/D(i)^2)
                 end
             elseif order == -3 %derivative order 1
                 if j == 1
                     output(i) = output(i) + F_coff(j)* D(i)^3/6
                 else
                     output(i) = output(i) - (F_coff(j)/D(i)^3)*sin((j-1)*D(i)) + (F_coff(j)/D(i)^2)*D(i)
                 end
             elseif order == -4 %derivative order 1
                 if j == 1
                     output(i) = output(i) + F_coff(j)* D(i)^4/24
                 else
                     output(i) = output(i) + (F_coff(j)/D(i)^4)*cos((j-1)*D(i)) - (F_coff(j)/D(i)^4) + F_coff(j)            
                 end
            end
        end
    end 
end

function output = test_der_1(F, D, F_coff)
    N = length(D);
    deg = leng(F_coff);
    output = zero(1,N);
    for i = 1:N
        for j = 1:deg
            output(i) = output(i) + F_coff(j)*cos(j*D(i));
        end
    end
end

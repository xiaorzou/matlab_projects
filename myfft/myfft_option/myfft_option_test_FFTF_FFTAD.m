%input
%   @M_A: degree of antiderivative  
%   @M_D: degree of derivative 
%{
M_A=6; 
M_D=6;
%}

function myfft_option_test_FFTF_FFTAD(M_A, M_D)
    addpath('D:/matlab/mylib')  
    addpath('D:/matlab/myfft/myfft_engine')  
    models={'Gauss','Merton','Kou','NIG','VG','TS'};
    modelTypes=[1,2,3,4,5,6];
    %modelTypes=[6];
    scenario = 1;
    for modelType = modelTypes 
        model = char(models(modelType))
        input_path=strcat('');
        input_fileName='input_test_FFTF_FFTAD_EXP.xlsx'; 
        [a, b] =xlsread(strcat(input_path,input_fileName),model);
        inputVar = containers.Map(b,a(scenario,:));
 
        [fftpower, left, right, ] = myfft_option_get_fft_para(scenario);
        N=2^fftpower;

        density_coef=myfft_density_coef(model, inputVar, left,right, N);
        %max(abs(F-F_new))
        
        FFTF = myfft_derivative(density_coef, M_D, left, right);
        file_fftf = ['output/fftf_' model '_power_' int2str(fftpower) '_right_' int2str(right) '.dat'];
        vnames  = {};
        for i = 1:M_D
            vnames = [vnames, strcat('deg',int2str(i))];
        end  
        mylib_writearray(vnames, transpose(FFTF), file_fftf)
        
        FFTAD = myfft_antiderivative(density_coef, M_A, left, right);
        file_fftad = ['output/fftad_' model '_power_' int2str(fftpower) '_right_' int2str(right) '.dat'];
        vnames  = {};
        for i = 1:M_A
            vnames = [vnames, strcat('deg',int2str(i))];
        end
        mylib_writearray(vnames, transpose(FFTAD), file_fftad)

        deg = 5;
        FFTExp = myfft_integral_exp(density_coef, deg, left, right);
        file_fftexp = ['output/fftexp_' model '_power_' int2str(fftpower) '_right_' int2str(right)  '.dat'];
        vnames  = {};
        for i = 1:deg
            vnames = [vnames, strcat('deg',int2str(i))];
        end
        mylib_writearray(vnames, transpose(FFTExp), file_fftexp);
    end
end

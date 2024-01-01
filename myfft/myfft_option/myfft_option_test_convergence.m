%input
%   @right:  [-right, right] effective range of q(x)
%   selected_power: 2^selected_power+1 points in plot/output 
%   samples of input:
%   right = 20;
%   selected_power = 6;
function myfft_option_test_convergence(right,selected_power)
    addpath('D:/matlab/mylib')  
    addpath('D:/matlab/myfft/myfft_engine')  
    models={'Gauss','Merton','Kou','NIG','VG','TS'};
    scenario = 1;
    %modelTypes=[1,2,3,4,5,6];
    modelTypes=[6];
    left = -right;
    power_base = 10;
    fftpowers = 10:20;
    selected_num = 2^selected_power;
    [x,CP_base] = myfft_get_grid(power_base, right);
    x = x((-selected_num:selected_num)+CP_base); %remove the last, so length(x)=2^power_base/2
    x = transpose(x);
    for modelType = modelTypes 
        model = char(models(modelType))
        input_path=strcat('');
        input_fileName='input_test_FFTF_FFTAD_EXP.xlsx'; 
        [a, b] =xlsread(strcat(input_path,input_fileName),model);
        inputVar = containers.Map(b,a(scenario,:));
        output = zeros(2*selected_num+1, length(fftpowers)+1);
        output_plot = zeros(2*selected_num+1, length(fftpowers)*2);
        counter = 2;
        output(:,1) = x;
        for fftpower = fftpowers
            N=2^fftpower;
            [D,CP] = myfft_get_grid(fftpower, right);
            pos_selected =  1:2^(fftpower-power_base):(N/2);
            plot_selected = (-selected_num:selected_num)+CP;
            output_plot(:, (counter-1)*2-1) = transpose(D(plot_selected));
            density_coef=myfft_density_coef(model, inputVar, left,right, N);
            FFTF = myfft_derivative(density_coef, 1, left, right);
            output_plot(:, (counter-1)*2) = transpose(FFTF(1,plot_selected));
            y = FFTF(1,pos_selected);
            y = y((-selected_num:selected_num)+CP_base);
            x_back = D(pos_selected);
            x_back = x_back((-selected_num:selected_num)+CP_base);
            output(:,counter)= transpose(y);
            diff = max(abs(transpose(x_back)-x))
            counter = counter + 1;
        end
        file_fftf = ['output/test_convergence_' model '_power_' int2str(fftpower) '_right_' int2str(right) '.dat'];
        file_plot = ['output/test_convergence_plot_' model '_power_' int2str(fftpower) '_right_' int2str(right) '.dat'];
        vnames  = {'x'};
        vnames_plot = {};
        for i = 1:length(fftpowers)
            vnames = [vnames, strcat('power_',int2str(fftpowers(i)))];
            vnames_plot = [vnames_plot, strcat('x_',int2str(fftpowers(i))), strcat('y_',int2str(fftpowers(i)))];
        end  
        mylib_writearray(vnames, output, file_fftf)
        mylib_writearray(vnames_plot, output_plot, file_plot)
    end
end

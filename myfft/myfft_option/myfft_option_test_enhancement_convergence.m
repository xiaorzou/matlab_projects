%input
%   @right:  [-right, right] effective range of q(x)
%   selected_power: 2^selected_power+1 points in plot/output 
%   samples of input:
%   right = 20;
%   selected_power = 6;
function myfft_option_test_enhancement_convergence()
    addpath('D:/matlab/mylib')  
    addpath('D:/matlab/myfft/myfft_engine')  
    addpath('D:/matlab/fourier_normalizer/fourier_normalizer_engine')  
    models={'Gauss','Merton','Kou','NIG','VG','TS'};
    scenario = 1;
    selected_power = 6;
    modelTypes=[1,2,3,4,5,6];
    %modelTypes=[6];
    right = 5;
    left = -right;
    b = 5;
    a = 4;
    c = 0.1;
    d = 0.1;
    n = 0;
    p = 8; %if p=10,then max(abs(h_coef(2^8,2^10)))<10^-12
    q = p+2;
    h_coef = fourier_normalizer_get_h_coef(a,b,c,d,n,p,q);
    h_value = fourier_normalizer_get_h_value(a,b,c,d,n,p,q,h_coef);
    x_prime = fourier_normalier_get_grid(2^p, 0, b);
    plot(x_prime,h_value)
    N_h = length(h_coef);
    Aones=ones(1,N_h);
    Aones(2:2:end)=-ones(1,N_h/2);
    h_coef_alt = h_coef.*Aones;
    N_h = 2*N_h;
    h_coef_new = zeros(1,N_h);
    h_coef_new(1:2:end) = h_coef_alt;
    clear h_coef Aones h_coef_alt
    power_base = 10;
    %fftpowers = 10:20;
    fftpowers = 10; % must not be less than p
    selected_num = 2^selected_power;
    [x,CP_base] = myfft_get_grid(power_base, right);
    x = x((-selected_num:selected_num)+CP_base); %remove the last, so length(x)=2^power_base/2
    x = transpose(x);
    for modelType = modelTypes 
        model = char(models(modelType))
        input_path=strcat('');
        input_fileName='input_test_FFTF_FFTAD_EXP.xlsx'; 
        [scen_value, scen_att] =xlsread(strcat(input_path,input_fileName),model);
        inputVar = containers.Map(scen_att,scen_value(scenario,:));
        output = zeros(2*selected_num+1, length(fftpowers)+1);
        output_enhance = zeros(2*selected_num+1, length(fftpowers)+1);
        output_plot = zeros(2*selected_num+1, length(fftpowers)*3);
        counter = 2;
        output(:,1) = x;
        output_enhance(:,1) = x;
        for fftpower = fftpowers    
            N=2^fftpower;
            [D,CP] = myfft_get_grid(fftpower, right);
            pos_selected =  1:2^(fftpower-power_base):(N/2);
            plot_selected = (-selected_num:selected_num)+CP;
            
            density_coef_Nplus=myfft_density_coef(model, inputVar, left,right, 2*N); %need double N to do convolution with h
            
            density_coef_N = density_coef_Nplus(1:N);
            FFTF = myfft_derivative(density_coef_N, 1, left, right);
            y = FFTF(1,pos_selected);
            y = y((-selected_num:selected_num)+CP_base);
            output(:,counter)= transpose(y);
            
            density_coef_normalized = myfft_coef_normalizer(density_coef_Nplus, h_coef_new, N);
            %{ 
            remove the following two lines to 'myfft_coef_normalizer'
            density_coef_normalized(1) = h_coef_new(1)*density_coef_N(1) + 0.5*dot(h_coef_new(2:N_h), density_coef_N(2:N_h));
            density_coef_normalized(2:N_h) = density_coef_normalized(2:N_h)+ 0.5*density_coef_N(1).*h_coef_new(2:N_h);
            %}
            
            FFTF_enhance = myfft_derivative(density_coef_normalized, 1, left, right);
            
            output_plot(:, (counter-1)*3-2) = transpose(D(plot_selected));
            output_plot(:, (counter-1)*3-1) = transpose(FFTF(1,plot_selected));
            output_plot(:, (counter-1)*3) = transpose(FFTF_enhance(1,plot_selected));

            y_enhance = FFTF_enhance(1,pos_selected);
            y_enhance = y_enhance((-selected_num:selected_num)+CP_base);
            output_enhance(:,counter)= transpose(y_enhance);

            counter = counter + 1;
            
            
        end
        file_fftf = ['output/test_enhancement_convergence_' model '_power_' int2str(fftpower) '_right_' int2str(right) '.dat'];
        file_plot = ['output/test_enhancement_convergence_plot_' model '_power_' int2str(fftpower) '_right_' int2str(right) '.dat'];
        vnames  = {'x'};
        vnames_plot = {};
        for i = 1:length(fftpowers)
            vnames = [vnames, strcat('power_',int2str(fftpowers(i)))];
            vnames_plot = [vnames_plot, strcat('x_',int2str(fftpowers(i))), strcat('y_',int2str(fftpowers(i))),strcat('y_enhance',int2str(fftpowers(i)))];
        end  
        mylib_writearray(vnames, output, file_fftf)
        mylib_writearray(vnames_plot, output_plot, file_plot)
    end
end

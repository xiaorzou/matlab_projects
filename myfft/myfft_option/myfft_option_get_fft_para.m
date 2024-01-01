
%input
%   @model: integer, example: scenario = 1
%output
%   @fftpower:  N=2^fftpower is used in fourier expansion with effective term N/2.
%   @Left: left end of effective range of density function.
%   @Right: right end of effective range of density function.
%   @M_A: up to M_A antiderivative of density function is handled. fourier expansion of CFD is
%   assoicated to FFTAD(1,:).
%   @M_D: up to M_D-1 derivative of density function is handled. fourier
%   expansion of density if assoicated to FFTF(1,:)

function  [fftpower, Left, Right] = myfft_option_get_fft_para(scenario)
    input_path=strcat('');
    input_myfft_fileName='input_myfft.xlsx'; 
    [values, index] =xlsread(strcat(input_path,input_myfft_fileName),'FFT'); 
    FFTVar = containers.Map(index,values(scenario,:));
    fftpower = FFTVar('fftpower');
    Left = FFTVar('Left');
    Right = FFTVar('Right'); 
end

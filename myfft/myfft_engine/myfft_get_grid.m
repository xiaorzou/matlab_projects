% 
% input:
%   @fftpower: 2^fftpower is the number of items in cos expansion of q(x)
%   @right:  [-right, right] defines the effective range of the density function q(x)
% output: 
%   @D: grid points, starting from left, length(D)=N/2+1
%   @CP: the center position,  associated to 0 value in grid. 


function [D,CP] = myfft_get_grid(fftpower, right)
    left = -right; %force symmetric, CP=N/4+1 is assoicated to the pos of 0
    N=2^fftpower;
    lambda=(right-left)/(N/2);
    D=(0:1:N/2)*lambda+left; 
    CP=N/4+1;
end


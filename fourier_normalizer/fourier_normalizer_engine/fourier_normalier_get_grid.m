% input
% left:  the left point of the  interval [left, right]
% right: the right point of the  interval [left, right]
% N:  the number of enevely divided interval grid points in the range (left-ll, right) 
% output: 
% the array of grid points to divide the interval (left-ll, right) with ll = right-left

function val = fourier_normalier_get_grid(N, left, right)
ll = right - left;
delta = 2*ll/N;
val = (0:1:(N-1))*delta+left-ll;
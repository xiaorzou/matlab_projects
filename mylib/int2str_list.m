function val = int2str_list( list_int )
%INT2STR_LIST Summary of this function goes here
%   localtion:  matlab/mylib
%   input:  
%       @list_int: a list of integers;   example [1,2,3]
%   output:
%       @val: assoicated list with char;   example ['1','2','3']
    val = strtrim(cellstr(num2str(list_int'))');
end


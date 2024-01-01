function mylib_writearray(vnames, input_array, filename_full)
    input_array = array2table(input_array, 'VariableNames',vnames);
    writetable(input_array,filename_full)    
end


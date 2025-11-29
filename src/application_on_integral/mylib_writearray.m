function mylib_writearray(vnames, input_array, filename_full)
    input_array = array2table(input_array, 'VariableNames',vnames);
    if exist(filename_full, 'file')
        delete(filename_full);
    end
    writetable(input_array,filename_full)
end


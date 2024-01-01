% the function to get h(x) in the doc.
function max_table = fourier_normalizer_test_get_h_w()
    fig_type = '.fig';
    flag_plot = true;
    p = 9;
    a = 2;
    ns = 0:5;
    test_type = 'full';
    if strcmp(test_type,'full')
        bs = [3,7,10,15];
        cs = [0.05,0.1,0.5,1];
        ds = [0.05,0.1,0.5,1];
        qs = [10,15,17];
    else
        bs = [3];
        cs = [0.1];
        ds = [0.1];
        qs = [10];
    end
    
    vnames = {'a', 'b', 'c', 'd','n','p','q', 'max_q'};
    max_record = zeros(length(ns)*length(qs)*length(bs)*length(cs)*length(ds),length(vnames));
    counter = 1;
    for n = ns
        for q = qs
            for bid = 1:length(bs)
                b=bs(bid);
                for cid = 1:length(cs)
                    c=cs(cid);
                    for did = 1:length(ds)
                        d = ds(did);
                        h_w_excel_file = ['output/h_a_b_c_d_n_p_q_' num2str(a) '_' num2str(b) '_' num2str(c) '_' num2str(d) '_' num2str(n) '_' num2str(p) '_' num2str(q)  '.xlsx'];
                        output_table = fourier_normalizer_get_h_w(a,b,c,d,n,p,q, flag_plot, fig_type);
                        writetable(output_table{1},h_w_excel_file)
                        max_record(counter,:)= [a,b,c,d,n,p,q,output_table{2}];
                        counter = counter + 1;
                    end
                end
            end
        end
    end
    h_w_excel_file = 'output/fourier_normalizer_test_get_h_w_output_max.xlsx';
    max_table = array2table(max_record, 'VariableNames',vnames);
    writetable(max_table,h_w_excel_file)
    
end






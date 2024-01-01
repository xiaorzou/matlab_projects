function fourier_normalizer_test_h_at_boundary()
    check_bdy_deg = 5;
    p = 9;
    a = 2;
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
    
    
    M = 2^(p);
    rows = length(bs)*length(cs)*length(ds)*length(qs)*length(check_bdy_deg);
    deravative_at_pi_deg_n = zeros(rows,check_bdy_deg+7);
    deravative_at_0_deg_n = zeros(rows,check_bdy_deg+7);
    counter = 1;
    Aones2=ones(1,M);
    Aones2(2:2:end)=-ones(1,M/2);
    for q = qs
        for n = 0:check_bdy_deg
            for bid = 1:length(bs)
                b=bs(bid);
                for cid = 1:length(cs)
                    c=cs(cid);
                    for did = 1:length(ds)
                        d = ds(did);
                        h_coef = fourier_normalizer_get_h_coef(a,b,c,d,n,p,q);
                        pib = pi*(0:1:M-1)/b;
                        pib = pib.*pib;
                        pib_2 = pib;
                        alternative  = -1;
                        deravative_at_pi_deg_n(counter,1:7) = [n,a,b,c,d,p,q];
                        deravative_at_0_deg_n(counter,1:7) = [n,a,b,c,d,p,q];
                        for i=1:check_bdy_deg
                            deravative_at_pi_deg_n(counter,7+i) = log10(abs(alternative*dot(h_coef,pib.*Aones2)));
                            deravative_at_0_deg_n(counter,7+i) = log10(abs(alternative*dot(h_coef,pib)));
                            pib = pib.*pib_2;
                            alternative = -alternative;
                        end
                        counter = counter + 1;
                    end
                end
            end
        end
    end
    vnames = {'n', 'a', 'b', 'c', 'd','p','q'};
    for deg = 1:check_bdy_deg
        vnames = [vnames , {strcat('deg',num2str(2*deg))}];
    end
    deravative_at_pi_table = array2table(deravative_at_pi_deg_n, 'VariableNames',vnames);
    deravative_at_0_table = array2table(deravative_at_0_deg_n, 'VariableNames',vnames);
    file_name_left = 'output/fourier_normalizer_test_h_at_boundary_left.dat';
    file_name_middle = 'output/fourier_normalizer_test_h_at_boundary_middle.dat';
    writetable(deravative_at_pi_table,file_name_left)
    writetable(deravative_at_0_table,file_name_middle)
    
    
end


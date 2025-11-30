
%{
the driver for the following paper
"An Application of the Trigonometric Interpolation on Second Order Linear ODE"
%}

function main_tiba_with_w()
    global cuf_off_para 
    
    figtype_fig = '.fig';
    plot_location = 'northwest';
    fourier_normalizer_const()
    s = 1;
    e = 3*s;
    flag_save_excel = true;
    flag_save_fig = false;
    smooth_PQR = true; %default
    cond_types = {'Neumann','Dirichlet', 'mix_1', 'mix_2'}; %in the paper
    qs = [6,7,8]; %in paper 
    degs = [2]; % in paper 
    etas = [0, 1, 2]; % in paper 
    paras =[0.5, 1, 2, 4]; %in paper
    sts = [e];  %%in the paper,  if include s,  will have warning message on Neumann type!    
    

    p_u = 0.1;
    p_v = 1.0;

    function val = fun_sol(x, deg, theta)    %only for x<=0 old:   change to x>=0
        val = x.^deg.*cos(theta*x);
    end

    function val = fun_der_sol(x, theta, deg)  
        if deg==1
            val = cos(theta*x)-theta*x.*sin(theta*x);
        else
            val = deg*x.^(deg-1).*cos(theta*x)-theta*x.^deg.*sin(theta*x);
        end 
    end

    function val = fun_der2_sol(x, deg, theta)  %only for x<=0
        if deg==1
            val = -2*theta*sin(theta*x)-theta^2*x.*cos(theta*x);
        else
            val = deg*(deg-1)*x.^(deg-2).*cos(theta*x)-2*deg*theta*x.^(deg-1).*sin(theta*x)-theta^2*x.^deg.*cos(theta*x);
        end
    end

    function val = fun_w(x, eta, st)
        mm = length(x);
        if eta==0
            val = ones(1, mm);
        else
            val = (x-st).^eta;
        end
    end

    function val = generate_r(x, eta, st, deg, theta, p_u, p_v)
        z =  fun_der2_sol(x, deg, theta);
        v = fun_sol(x, deg, theta);  
        u = fun_der_sol(x, theta, deg);
        w = fun_w(x,eta, st);
        val = w.*z-p_v*v-p_u*u; 
    end


    function val = fun_q(X, p_v)
        val = p_v*ones(1,length(X));
    end

    function val = fun_p(X,p_u)
        val = p_u*ones(1,length(X));
    end


    output_max_error = zeros(length(cond_types)*length(qs)*length(sts)*length(paras)*length(degs)*length(etas),7);
    counter = 1;
    for q = qs
        for st = sts
            for eta = etas
                for para = paras
                    for deg = degs
                        p = q-1; %p should be smaller than q,  
                        theta = para*pi;

                        M = 2^q;
                        N = 2*M;
                        n = 2^p; % we should make n as large as possible,  as such, p=q-1 should be always true
                        lambda = (e-s)/n;
                        m = (M-n)/2;
                        delta = lambda*m;
                        o = s - delta;
                        b = e-s + 2*delta;
                        x_N = -b + o + lambda*(0:N-1);
                        x_R = x_N(M+1:end);
                        co_R = fourier_normalizer_cut_off(x_R, s-delta,s,e,e+delta,cuf_off_para);
                        z_1_true = fun_sol(s, deg, theta);
                        z_2_true = fun_der_sol(s, theta, deg);
                        z_3_true = fun_sol(e, deg, theta);
                        z_4_true = fun_der_sol(e, theta, deg);
                        v_s =  fun_sol(s, deg, theta);

                        for jj = 1:length(cond_types)
                            cond_type = cond_types(jj);

                            if strcmp(cond_type, 'Dirichlet')
                                A1 = [1,0,0,0];
                                A2 = [0,0,1,0]; 
                            elseif strcmp(cond_type, 'Neumann')
                                A1 = [1,0,0,0];
                                A2 = [0,1,0,0]; 
                            elseif strcmp(cond_type, 'mix_1')
                                A1 = [1,0,0,0];
                                A2 = [0,0,0,1];
                            elseif strcmp(cond_type, 'mix_2')
                                A1 = [1,1,0,0];
                                A2 = [0,0,1,1];
                            end
                            alpha = A1(1)*z_1_true+A1(2)*z_2_true+A1(3)*z_3_true+A1(4)*z_4_true;
                            beta = A2(1)*z_1_true+A2(2)*z_2_true+A2(3)*z_3_true+A2(4)*z_4_true;
                            Q = fun_q(x_R, p_v);
                            P = fun_p(x_R,p_u);                            
                            R = generate_r(x_R,eta, st, deg, theta, p_u, p_v);
                            
                            if smooth_PQR
                                Q = Q.*co_R;
                                P = P.*co_R;
                                R = R.*co_R;
                            end
                                
                            W = fun_w(x_R,eta, st);
                            
                            V_linear =  engine_tiba_with_w(b, W, R, Q, P, alpha, beta, m, n, cond_type,  A1,A2, v_s);   

                            x_R_plusb = zeros(1, M+1);
                            x_R_plusb(1:M) = x_R;
                            x_R_plusb(M+1) = b;
                            v_true_this = fun_sol(x_R_plusb, deg, theta);
                            pos_this_ = find(x_R_plusb>=s & x_R_plusb<=e);                                                                                    
                            max_error_linear = max(abs(v_true_this(pos_this_)'-V_linear(pos_this_)));                                                                                          
                            output_max_error(counter,:) = [q, st, eta, theta, deg, jj, max_error_linear];
                            counter = counter +1;
                            if flag_save_fig && q == 7 && strcmp(cond_type, 'Dirichlet')
                                xlabel_this = '$x$';  
                                ylabel_this = '$y$';
                                filename_fig = strcat('output/tiba_linear_enh', '_q_', num2str(q), '_eta_', num2str(eta), '_deg_', num2str(deg), '_theta_', num2str(round(theta)), '_st_', num2str(st), '_smooth_PWR_', num2str(smooth_PQR), figtype_fig);
                                title = strcat('$y_{tiba}$ vs $y_{true}$ over $[0,b]$');
                                legend_y = '$y_{tiba}$';
                                legend_z = '$y_{true}$';                                
                                plot_latex_2(filename_fig,x_R_plusb, V_linear',  v_true_this, title, legend_y, legend_z,  xlabel_this, ylabel_this, plot_location);                    
                            end
                        end
                    end
                end
            end
        end
    end
    if flag_save_excel
        fourier_normalizer_app_ode_deg_2_max_error_file = ['output/main_tiba_with_w_smooth_PQR_', num2str(smooth_PQR),  '.xlsx'];
        vnames = {'q', 'singular','eta', 'theta', 'deg','cond_type','linear_err'};
        mylib_writearray(vnames, output_max_error, fourier_normalizer_app_ode_deg_2_max_error_file);
    end
end

    

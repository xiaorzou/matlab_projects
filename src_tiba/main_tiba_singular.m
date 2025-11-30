
%{
the driver for the following paper
"Trigonometric Interpolation Based Based Approach for Second Order ODE with
Mixed Boundary Conditions"
%}

function main_tiba_singular()
    global cuf_off_para 
    
    figtype_fig = '.fig';
    plot_location = 'northwest';
    fourier_normalizer_const()
    s = 1;
    e = 3*s;
    flag_save_excel = true;
    flag_save_fig = true;
    smooth_PQR = true;
    %cond_types = {'Neumann','Dirichlet', 'mix_1', 'mix_2'}; 
    cond_types = {'Neumann','Dirichlet'}; % in the paper
    p_u = 1.0;
    p_v = 0.0;

    function val = fun_sol(x, theta)    %only for x<=0 old:   change to x>=0
        val = cos(theta*x);
    end

    function val = fun_der_sol(x, theta)  
       val = -theta*sin(theta*x);
    end

    function val = fun_der2_sol(x, theta)  %only for x<=0
        val = -theta^2*cos(theta*x);
    end

    function val = fun_w(x, st)
        val = (x-st);
    end

    function val = generate_r(x, st, theta, p_u, p_v)
        z =  fun_der2_sol(x, theta);
        v = fun_sol(x, theta);  
        u = fun_der_sol(x, theta);
        w = fun_w(x, st);
        %val = z-p_vv*v.^2-p_v*v-p_u*u;
        val = w.*z-p_v*v-p_u*u; %only linear term!
    end


    function val = fun_q(X, p_v)
        %val = p_v*ones(1,length(X)).*co_R;
        val = p_v*ones(1,length(X));
    end

    function val = fun_p(X,p_u)
        %val = p_u*ones(1,length(X)).*co_R;
        val = p_u*ones(1,length(X));
    end

    qs = [6,7,8]; %default 

    paras = [1,2,3]; %default 
    %sts = [s,e];
    st = e;
    output_max_error = zeros(length(qs)*length(paras),5);
    counter = 1;
    for q = qs
        for para = paras            
            %q = 7;
            p = q-1; %p should be smaller than q,  
            %theta = (0.5+ para)*pi;
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
            z_1_true = fun_sol(s,  theta);
            z_2_true = fun_der_sol(s, theta);
            z_3_true = fun_sol(e, theta);
            z_4_true = fun_der_sol(e, theta);
            v_s =  fun_sol(s, theta);
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
                %A1 = [1,0,0,0];
                %A2 = [0,0,1,0]; 

                alpha = A1(1)*z_1_true+A1(2)*z_2_true+A1(3)*z_3_true+A1(4)*z_4_true;
                beta = A2(1)*z_1_true+A2(2)*z_2_true+A2(3)*z_3_true+A2(4)*z_4_true;
                Q = fun_q(x_R, p_v);

                P = fun_p(x_R,p_u);

                R = generate_r(x_R, st, theta, p_u, p_v);

                if smooth_PQR
                    Q = Q.*co_R;
                    P = P.*co_R;
                    R = R.*co_R;
                end

                W = fun_w(x_R, st);                            

                V_linear =  engine_tiba_with_w(b, W, R, Q, P, alpha, beta, m, n, cond_type,  A1,A2, v_s);   
                x_R_plusb = zeros(1, M+1);
                x_R_plusb(1:M) = x_R;
                x_R_plusb(M+1) = b;
                v_true_this = fun_sol(x_R_plusb, theta);
                pos_this_ = find(x_R_plusb>=s & x_R_plusb<e);
                max_error_linear = max(abs(v_true_this(pos_this_)'-V_linear(pos_this_)));                   
                output_max_error(counter,:) = [q, st, theta, jj, max_error_linear];
                counter = counter +1;
                if flag_save_fig && q==7 &&  strcmp(cond_type, 'Dirichlet')
                    xlabel_this = '$x$';  
                    ylabel_this = '$y$';
                    filename_fig = strcat('output/main_tiba_signular', '_q_', num2str(q), '_theta_', num2str(round(theta)), '_', cond_type{1} ,  figtype_fig);
                    title = strcat('$y_{tiba}$ vs $y_{true}$ over $[0,b]$');
                    legend_y = '$y_{tiba}$';
                    legend_z = '$y_{true}$';                                
                    plot_latex_2(filename_fig,x_R_plusb(pos_this_), V_linear(pos_this_)'./W(pos_this_),  v_true_this(pos_this_)./W(pos_this_), title, legend_y, legend_z,  xlabel_this, ylabel_this, plot_location);                    
                end
            end
        end
    end
    if flag_save_excel
        fourier_normalizer_app_ode_deg_2_max_error_file = ['output/main_tiba_signular_smooth_PQR_', num2str(smooth_PQR), '.xlsx'];
        vnames = {'q', 'singular', 'theta', 'bdy' , 'linear_err'};
        mylib_writearray(vnames, output_max_error, fourier_normalizer_app_ode_deg_2_max_error_file);
    end
end

    

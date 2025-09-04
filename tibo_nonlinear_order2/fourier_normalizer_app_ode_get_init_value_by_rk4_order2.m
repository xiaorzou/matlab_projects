
%{
fun:  f(x,y,y')
[s,e]:  boundary point
two combination of init/boundary conditions:
alpha = A(1,1)*z_1_true+A(1,2)*z_2_true+A(1,3)*z_3_true+A(1,4)*z_4_true;
beta = A(1,2)*z_1_true+A(2,2)*z_2_true+A(2,3)*z_3_true+A(2,4)*z_4_true;
inits:  inital guess for (y(s),y'(s))
zopt: opt value for (y(s),y'(s))
myobjopt: opt error
%}


%function [zopt, myobjopt, X_t] = fourier_normalizer_app_ode_get_init_value_by_rk4_order2(fun, s, e,  alpha, beta, A1, A2, J, inits, delta)
%function [zopt, myobjopt, X_t] = fourier_normalizer_app_ode_get_init_value_by_rk4_order2(fun, s, e,  greeks, AA, J, inits, delta, maxiter, constrainttolerance)
function [zopt, myobjopt, X_t] = fourier_normalizer_app_ode_get_init_value_by_rk4_order2(fun, s, e,  greeks, AA, J, inits, maxiter, constrainttolerance)
    addpath('D:/matlab/fourier_normalizer/fourier_normalizer_engine')  
    addpath('D:/matlab/cos_approx/cos_approx_engine')  
    %addpath('D:/matlab/mylib') 
    lambda = (e-s)/J;
    delta_shock = 0.001; %this delta is different from TIBO. 
    t_vector = s + (0:J)*lambda;
    pos = 1;
    alpha = greeks(1);
    beta = greeks(2);
    flag_fix_1 = false;
    flag_fix_2 = false;
    if AA(1,3)^2+AA(1,4)^2==0 && AA(2,3)^2+AA(2,4)^2==0
        AAA = [AA(1,1:2); AA(2,1:2)];
        zopt = AAA\[alpha;beta];
        myobjopt = 0;
        [X_t, ~] = runge_kutta_ode_order_k(fun, t_vector,zopt, pos, []);
        return
    end
    if AA(1,2)^2+AA(1,3)^2+AA(1,4)^2==0
        flag_fix_1 = true;
        inits(1)=alpha/AA(1,1);
    end
    if AA(2,1)^2+AA(2,3)^2+AA(2,4)^2==0
        flag_fix_2 = true;
        inits(2)=beta/AA(2,2);
    end
  
    %{
    t_11 = A(1,1) + A(1,2);
    t_12 = A(1,1)*s + A(1,2)*e + A(1,3)+A(1,4);
    t_21 = A(2,1) + A(2,2);
    t_22 = A(2,1)*s + A(2,2)*e + A(2,3)+A(2,4);
    %}
    
    flag_display = false;
    %flag_display = true;
    %maxiter = 10000;
    %constrainttolerance = 10^-10;
    
    if flag_fix_1
        z_init =  inits(2);
    elseif flag_fix_2
        z_init =  inits(1);
    else
        z_init = inits;
    end
   
    zlast = [];
    myobjlast = [];
    glast = [];
    hlast = [];
    dmyobjlast = [];
    dglast = [];
    dhlast = [];
    lb = [];
    ub = [];
    X_t = [];
    function [myobj,g,h,dmyobj,dg,dh] = objcon(z_M) 
        
        if flag_fix_1
            z_1 = inits(1);
            z_2 = z_M(1);
            z_2_p1 = z_2+delta_shock;
            z_2_m1 = z_2-delta_shock;
            [X_t, ~] = runge_kutta_ode_order_k(fun, t_vector,[z_1, z_2], pos, []);
            [X_2p, ~] = runge_kutta_ode_order_k(fun, t_vector,[z_1, z_2_p1], pos, []);
            [X_2m, ~] = runge_kutta_ode_order_k(fun, t_vector,[z_1, z_2_m1], pos, []);
            z3_2 = (X_2p(1, end) - X_2m(1, end))/(2*delta_shock);
            z4_2 = (X_2p(2, end) - X_2m(2, end))/(2*delta_shock);
            z3_1 = 0;
            z4_1 = 0;
            z_3 = X_t(1,end);
            z_4 = X_t(2,end);
        elseif flag_fix_2
            z_1 = z_M(1);
            z_2 = inits(2);
            z_1_p1 = z_1+delta_shock;
            z_1_m1 = z_1-delta_shock;
            [X_t, ~] = runge_kutta_ode_order_k(fun, t_vector,[z_1, z_2], pos, []);
            [X_1p,~] = runge_kutta_ode_order_k(fun, t_vector,[z_1_p1, z_2], pos, []);
            [X_1m, ~] = runge_kutta_ode_order_k(fun, t_vector,[z_1_m1, z_2], pos, []);
            z3_1 = (X_1p(1, end) - X_1m(1, end))/(2*delta_shock);
            z4_1 = (X_1p(2, end) - X_1m(2, end))/(2*delta_shock);
            z3_2 = 0;
            z4_2 = 0;
            z_3 = X_t(1,end);
            z_4 = X_t(2,end);

        else
            z_1 = z_M(1);
            z_2 = z_M(2);
            z_1_p1 = z_1+delta_shock;
            z_1_m1 = z_1-delta_shock;
            z_2_p1 = z_2+delta_shock;
            z_2_m1 = z_2-delta_shock;
            [X_t, ~] = runge_kutta_ode_order_k(fun, t_vector,[z_1, z_2], pos, []);
            [X_1p, ~] = runge_kutta_ode_order_k(fun, t_vector,[z_1_p1, z_2], pos, []);
            [X_1m, ~] = runge_kutta_ode_order_k(fun, t_vector,[z_1_m1, z_2], pos, []);
            [X_2p, ~] = runge_kutta_ode_order_k(fun, t_vector,[z_1, z_2_p1], pos, []);
            [X_2m, ~] = runge_kutta_ode_order_k(fun, t_vector,[z_1, z_2_m1], pos, []);
            z3_1 = (X_1p(1, end) - X_1m(1, end))/(2*delta_shock);
            z3_2 = (X_2p(1, end) - X_2m(1, end))/(2*delta_shock);
            z4_1 = (X_1p(2, end) - X_1m(2, end))/(2*delta_shock);
            z4_2 = (X_2p(2, end) - X_2m(2, end))/(2*delta_shock);
            z_3 = X_t(1,end);
            z_4 = X_t(2,end);
        end
        part1 = -(alpha- AA(1,1)*z_1-AA(1,2)*z_2-AA(1,3)*z_3-AA(1,4)*z_4);
        part2 = -(beta- AA(2,1)*z_1-AA(2,2)*z_2-AA(2,3)*z_3-AA(2,4)*z_4);
        f1 = part1*(AA(1,1)+AA(1,3)*z3_1+AA(1,4)*z4_1)+ part2*(AA(2,1)+AA(2,3)*z3_1+AA(2,4)*z4_1);
        f2= part1*(AA(2,1)+AA(2,3)*z3_2+AA(2,4)*z4_2)+ part2*(AA(2,1)+AA(2,3)*z3_2+AA(2,4)*z4_2);
        myobj = 0.5*part1^2 + 0.5*part2^2;
        myobj = myobj';
        if flag_fix_1
            dmyobj = f2;
        elseif flag_fix_2
            dmyobj = f1;
        else
            dmyobj = [f1;f2];
        end
        g = [];
        h = [];
        dg = [];
        dh = [];
    end
    if flag_display
        options = optimoptions(@fmincon, ...
            'Algorithm', 'interior-point',...
            'Display', 'Iter-detailed',...
            'GradObj','on',...
            'MaxIter', maxiter,...
            'MaxFunEvals',maxiter,...
            'TolFun', constrainttolerance,...
            'GradConstr','on');
    else
        options = optimoptions(@fmincon, ...
            'Algorithm', 'interior-point',...
            'GradObj','on',...
            'MaxIter', maxiter,...
            'MaxFunEvals',maxiter,...
            'TolFun', constrainttolerance,...
            'GradConstr','on');
    end
        

    function [myobj,dmyobj] = obj(z)
        if ~isequal(z,zlast)
            [myobjlast, glast, hlast, dmyobjlast, dglast, dhlast] = objcon(z);
            zlast = z;
        end
        myobj = myobjlast;
        dmyobj = dmyobjlast;
    end
    function [g,h,dg,dh] = con(z)
        if ~isequal(z,zlast)
            [myobjlast, glast, hlast, dmyobjlast, dglast, dhlast] = objcon(z);
            zlast = z;
        end
        g = glast;
        h = hlast;
        dg = dglast;
        dh = dhlast; 
    end
    [zopt, myobjopt,iter, opt] = fmincon(@obj, z_init, [],[], [], [], lb, ub, @con, options);
    if flag_fix_1
        zopt = [inits(1), zopt];
    elseif flag_fix_2
        zopt = [zopt, inits(2)];
    end
end

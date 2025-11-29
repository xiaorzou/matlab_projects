
%{
 https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods
%}

function y = runge_kutta_ode_orde_1(fun, x_vector,init, pos, y_true)
    N = length(x_vector);
    y = zeros(1, N);
    y(pos) = init;

    flag_use_true = true;
    if isempty(y_true)
        flag_use_true = false;
    end        
    for i = (pos+1):N
        h = x_vector(i) - x_vector(i-1);
        if flag_use_true
            y(i) = runge_kutta_ode_orde_1_step(fun, x_vector(i-1), y_true(i-1), h);
        else    
            y(i) = runge_kutta_ode_orde_1_step(fun, x_vector(i-1), y(i-1), h);
        end
    end
    for i = (pos-1):-1:1
        h = x_vector(i) - x_vector(i+1);
        if flag_use_true
            y(i) = runge_kutta_ode_orde_1_step(fun, x_vector(i+1), y_true(i+1), h);
        else
            y(i) = runge_kutta_ode_orde_1_step(fun, x_vector(i+1), y(i+1), h);
        end
    end
end

function y_next = runge_kutta_ode_orde_1_step(fun, x_n,y_n,h)
    k_1 = fun([x_n; y_n]);
    k_2 = fun([x_n+h/2; y_n+k_1*h/2]);
    k_3 = fun([x_n+h/2; y_n+k_2*h/2]);
    k_4 = fun([x_n+h; y_n+k_3*h]);
    y_next = y_n+ h/6*(k_1+2*k_2+2*k_3+k_4);
end
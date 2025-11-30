
%{
 https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods

y^(k)(t) = fun(t;y^(k-1), y^(k-2), ..., y)

y=x1
y'=x2
y''=x3
...
y^(k-2) = x(k-1)
y^(k-1) = x(k)

x1' = x2
x2' = x3
x(k-1)'= x(k)
x(k)'= fun(t;x(k),...,x(1))

one can derive
X' = F(X;t):=(x2,x3,...,x(k),fun(t;x(1),..., x(k))

inits = (x1(t(pos)), x2(t(pos)),...,x(k-1)(t(pos)),xk(t(pos)))
pos:  start point position
t_vector:  grid nodes
X_true over t_vector

y^(k)(x) =fun_grid(x,y', ...,y^(k-1))

%}

function [X, z] = runge_kutta_ode_order_k(fun_grid, t_vector,inits, pos, X_true)
    N = length(t_vector);
    dim = length(inits);
    X = zeros(N, dim);
    X(pos, :) = inits;
    flag_use_true = true;
    if isempty(X_true)
        flag_use_true = false;
    end        
    function val = fun_vector(t_n, X_n)
        val = [X_n(2:end), fun_grid([t_n, X_n]')];
    end
            
    for i = (pos+1):N
        h = t_vector(i) - t_vector(i-1);
        if flag_use_true
            X(i,:) = runge_kutta_ode_orde_2_step(@fun_vector, t_vector(i-1), X_true(i-1,:),h);
        else    
            X(i,:) = runge_kutta_ode_orde_2_step(@fun_vector, t_vector(i-1), X(i-1,:), h);           
        end
    end
    for i = (pos-1):-1:1
        h = t_vector(i) - t_vector(i+1);
        if flag_use_true
            X(i,:) = runge_kutta_ode_orde_2_step(@fun_vector, t_vector(i+1), X_true(i+1,:), h);
        else
            X(i,:) = runge_kutta_ode_orde_2_step(@fun_vector, t_vector(i+1), X(i+1,:), h);
        end
    end
    X = X';
    X_t = [t_vector; X];
    z = fun_grid(X_t);
end

function x_next = runge_kutta_ode_orde_2_step(fun_vector, t_n, X_n, h)
    k_1 = fun_vector(t_n, X_n);
    k_2 = fun_vector(t_n+h/2, X_n+h*k_1/2);
    k_3 = fun_vector(t_n+h/2, X_n+h*k_2/2);
    k_4 = fun_vector(t_n+h, X_n+h*k_3);
    x_next = X_n+ h/6*(k_1+2*k_2+2*k_3+k_4);
end

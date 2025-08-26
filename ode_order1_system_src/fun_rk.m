function val = fun_rk(x, U, x_2N_R, co_2N_R, co_2N_der_R, test_type, theta,deg_x, ps, qs, c3)
    XU = [x,U];
    val = fun4rk(XU, x_2N_R, co_2N_R, co_2N_der_R, test_type, theta,deg_x, ps, qs, c3);
end


function string_simulation_template01()
    num_masses = 10;
    total_mass = 1;
    tension_force = 1;
    string_length = 1;
    damping_coeff = 0.01;
    dx = string_length/(num_masses+1);
    amplitude_Uf = 0.1;
    omega_Uf = 10;
    %list of x points (including the two endpoints)
    xlist = linspace(0,string_length,num_masses+2);
    Uf_func = @(t_in) amplitude_Uf*cos(omega_Uf*t_in);
    dUfdt_func = @(t_in) -omega_Uf*amplitude_Uf*sin(omega_Uf*t_in);
    %generate the struct
    string_params = struct();
    string_params.n = num_masses;
    string_params.M = total_mass;
    string_params.Uf_func = Uf_func;
    string_params.dUfdt_func = dUfdt_func;
    string_params.Tf = tension_force;
    string_params.L = string_length;
    string_params.c = damping_coeff;
    string_params.dx = dx;
    %load string_params into rate function
    my_rate_func = @(t_in,V_in) string_rate_func01(t_in,V_in,string_params);
    %initial conditions
    U0 = 0;
    dUdt0 = 0;
    V0 = [U0*ones(num_masses,1);dUdt0*ones(num_masses,1)];
    tspan = [0,10];
    %run the integration
    [tlist,Vlist] = explicit_RK_variable_step_integration(my_rate_func,tspan,V0, 0.01, rk_method("dormandprince"), 4, 10^-8);
    %your code to generate an animation of the system
    plot_system(xlist, tlist, Vlist, Uf_func)
end
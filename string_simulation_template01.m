function string_simulation_template01()

    %define location and filename where video will be stored
    input_fname = "modes_and_stuff.avi";
    %create a videowriter, which will write frames to the animation file
    writerObj = VideoWriter(input_fname);
    %must call open before writing any frames
    open(writerObj);

    for sweep_masses = 2:4:18
        for sweep_mode = 1:ceil(sweep_masses/3)
    num_masses = sweep_masses;
    total_mass = 1;
    tension_force = 1;
    string_length = 1;
    damping_coeff = 0.01;
    dx = string_length/(num_masses+1);
    %generate the struct
    string_params = struct();
    string_params.n = num_masses;
    string_params.M = total_mass;
    string_params.Tf = tension_force;
    string_params.L = string_length;
    string_params.c = damping_coeff;
    string_params.dx = dx;
    
    % % Calculate the Mode Shapes and Resonate Frequencies
    % omega = damping_coeff*((pi*n_int))/string_length;
    % Xn_x = Bn * sind(omega*x);
    w = sweep_mode;

    % Calculate Modal Analysis

    [M_mat, K_mat, ~] = construct_2nd_order_matrices(string_params);
    [Ur_mat,lambda_mat] = eig(K_mat,M_mat);
    omegas = sqrt(lambda_mat);
    amplitude_Uf = 0.05;
    omega_Uf = omegas(w,w);
    %list of x points (including the two endpoints)
    xlist = linspace(0,string_length,num_masses+2);
    Uf_func = @(t_in) amplitude_Uf*sin(omega_Uf*t_in);
    dUfdt_func = @(t_in) omega_Uf*amplitude_Uf*cos(omega_Uf*t_in);
    string_params.Uf_func = Uf_func;
    string_params.dUfdt_func = dUfdt_func;

    % Run Solver

    %load string_params into rate function
    my_rate_func = @(t_in,V_in) string_rate_func01(t_in,V_in,string_params);
    %initial conditions
    U0 = 0;
    dUdt0 = 0;
    V0 = [U0*ones(num_masses,1);dUdt0*ones(num_masses,1)];
    tspan = [0,((sweep_mode+1)*2*2*pi)/omega_Uf];

    %run the integration
    [tlist,Vlist] = explicit_RK_variable_step_integration(my_rate_func,tspan,V0, 0.01, rk_method("fehlberg"), 4, 10^-5);
    %your code to generate an animation of the system
    Vpred = (Ur_mat(:,w)*cos(omega_Uf*tlist))';
    plot_system(xlist, tlist, Vlist, Uf_func, Vpred, string_params, "# Masses = "+string(string_params.n)+"; Mode # = "+string(w), writerObj)
        end
    end
    close(writerObj)
end
function string_simulation_template01()
    num_masses = 10;
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
    string_params.wave_speed = sqrt(tension_force/(total_mass/string_length));
    
    % Predicted
    n = 2;
    freq = string_params.wave_speed*pi*n/string_params.L;
    continuous_x = linspace(0,string_length,1000);
    shape = @(x) sin(pi*n/string_params.L*x);
    figure;
    plot(continuous_x, shape(continuous_x), DisplayName="analytical")
    hold on;

    % Calculate Modal Analysis
    for num_masses = 2:2:20
        string_params.n = num_masses;
        [M_mat, K_mat, ~] = construct_2nd_order_matrices(string_params);
        [Ur_mat,lambda_mat] = eig(K_mat,M_mat);
        omegas = sqrt(lambda_mat);
        omega_Uf = omegas(n,n);
        Vpred = (Ur_mat(:,n)*cos(omega_Uf))';
        xlist = linspace(0,string_length,num_masses+2);
        mode_shape = [0,Vpred,0];
        V = mode_shape*(norm(shape(xlist))/norm(mode_shape));
        s = sign(mode_shape(2))*sign(shape(xlist(2)));
        plot(xlist, s*V, DisplayName=string(num_masses))
    end
    legend()
end
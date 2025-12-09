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
    n = 4;
    freq = string_params.wave_speed*pi*n/string_params.L;
    omega_freq = pi*n/string_params.L;
    continuous_x = linspace(0,string_length,1000);
    shape = @(x) sin(pi*n/string_params.L*x);
    figure;
    plot(continuous_x, shape(continuous_x), DisplayName="analytical")
    hold on;

    masses = 4:2:40;
    frequencies = zeros(40,length(masses));
    % Calculate Modal Analysis
    for i = 1:length(masses)
        num_masses = masses(i);
        string_params.n = num_masses;
        dx = string_length/(num_masses+1);
        string_params.dx = dx;

        [M_mat, K_mat, ~] = construct_2nd_order_matrices(string_params);
        [Ur_mat,lambda_mat] = eig(K_mat,M_mat);
        omegas = sqrt(lambda_mat);
        omega_Uf = omegas(n,n);
        d = diag(omegas);
        frequencies(1:length(d),i) = d;%/(2*pi);
        Vpred = (Ur_mat(:,n)*cos(omega_Uf))';
        xlist = linspace(0,string_length,num_masses+2);
        mode_shape = [0,Vpred,0];
        V = mode_shape*(norm(shape(xlist))/norm(mode_shape));
        s = sign(mode_shape(2))*sign(shape(xlist(2)));
        plot(xlist, s*V, DisplayName=string(num_masses))
    end
    legend()

    figure;
    semilogy(masses,abs(frequencies(n,:)-omega_freq)/omega_freq,".")
    xlabel("Masses")
    ylabel("Error Predicted Frequency (%)")

    omega_freqs = pi*(1:40)/string_params.L;
    figure;
    plot(1:40, omega_freqs, "." ,DisplayName="actual")
    hold on
    plot(1:40, frequencies(:,masses==4), "." ,DisplayName="4 masses");
    plot(1:40, frequencies(:,masses==8), "." ,DisplayName="8 masses");
    plot(1:40, frequencies(:,masses==16), "." ,DisplayName="16 masses");
    plot(1:40, frequencies(:,masses==32), "." ,DisplayName="32 masses");
    xlabel("# Harmonic")
    ylabel("Frequency (rad/s)")
    legend();
end
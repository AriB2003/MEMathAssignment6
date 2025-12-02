function string_simulation_template02()
    num_masses = 200;
    total_mass = 1;
    tension_force = 1;
    string_length = 1;
    damping_coeff = 0;
    dx = string_length/(num_masses+1);
    %generate the struct
    string_params = struct();
    string_params.n = num_masses;
    string_params.M = total_mass;
    string_params.Tf = tension_force;
    string_params.L = string_length;
    string_params.c = damping_coeff;
    string_params.dx = dx;
    string_params.pulse_width = 0.2;
    string_params.wave_speed = sqrt(tension_force/(total_mass/string_length));
    
    % Calculate Modal Analysis

    [M_mat, K_mat, ~] = construct_2nd_order_matrices(string_params);
    [Ur_mat,lambda_mat] = eig(K_mat,M_mat);
    omegas = sqrt(lambda_mat);

    w = 1;
    amplitude_Uf = 0.05;
    omega_Uf = omegas(w,w);
    %list of x points (including the two endpoints)
    xlist = linspace(0,string_length,num_masses+2);
    w = string_params.pulse_width;
    h = 0.2;
    Uf_func = @(t_in) b_spline_pulse(t_in, w,h);
    dUfdt_func = @(t_in) b_spline_pulse_derivative(t_in,w,h);
    string_params.Uf_func = Uf_func;
    string_params.dUfdt_func = dUfdt_func;

    % Run Solver

    %load string_params into rate function
    my_rate_func = @(t_in,V_in) string_rate_func01(t_in,V_in,string_params);
    %initial conditions
    U0 = 0;
    dUdt0 = 0;
    V0 = [U0*ones(num_masses,1);dUdt0*ones(num_masses,1)];
    tspan = [0,5];

    %run the integration
    [tlist,Vlist] = explicit_RK_variable_step_integration(my_rate_func,tspan,V0, 0.01, rk_method("fehlberg"), 4, 10^-6);
    %your code to generate an animation of the system
    % Vpred = (Ur_mat(:,w)*cos(omega_Uf*tlist))';
    plot_system(xlist, tlist, Vlist, Uf_func, 0, string_params)
end

%triangle pulse function
%INPUTS:
%t: current time
%w: width of pulse (starts at t=0, ends at t=h)
%h: height of pulse
%OUTPUTS:
%res: pulse evaluated at t
function res = triangle_pulse(t,w,h)
t = t*(2/w);
res = 1-min(1*abs(t-1),1);
res = h*res;
end
%triangle pulse function (derivative)
%INPUTS:
%t: current time
%w: width of pulse (starts at t=0, ends at t=h)
%h: height of pulse
%OUTPUTS:
%res: derivative of pulse evaluated at t
function res = triangle_pulse_derivative(t,w,h)
t = t*(2/w);
res = -sign(t-1).*(abs(t-1)<1);
res = (2*h/w)*res;
end
%b-spline pulse function
%INPUTS:
%t: current time
%w: width of pulse (starts at t=0, ends at t=h)
%h: height of pulse
%OUTPUTS:
%res: pulse evaluated at t
function res = b_spline_pulse(t,w,h)
t = 4*t/w;
b3 = (0<=t).*(t<1).*(t.^3)/4;
t = t-1;
b2 = (0<=t).*(t<1).*(-3*t.^3+3*t.^2+3*t+1)/4;
t = t-1;
b1 = (0<=t).*(t<1).*(3*t.^3-6*t.^2+4)/4;
t = t-1;
b0 = (0<=t).*(t<1).*(-t.^3+3*t.^2-3*t+1)/4;
res = h*(b0+b1+b2+b3);
end
%b-spline pulse function (derivative)
%INPUTS:
%t: current time
%w: width of pulse (starts at t=0, ends at t=h)
%h: height of pulse
%OUTPUTS:
%res: derivative of pulse evaluated at t
function res = b_spline_pulse_derivative(t,w,h)
t = 4*t/w;
b3 = (0<=t).*(t<1).*(3*t.^2)/4;
t = t-1;
b2 = (0<=t).*(t<1).*(-9*t.^2+6*t+3)/4;
t = t-1;
b1 = (0<=t).*(t<1).*(9*t.^2-12*t)/4;
t = t-1;
b0 = (0<=t).*(t<1).*(-3*t.^2+6*t-3)/4;
res = (4*h/w)*(b0+b1+b2+b3);
end
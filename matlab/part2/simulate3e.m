%% Clear
clear; close all;
to_plot = true;
to_print = false;
%% Constants
% Global constants 
V_a = 580 / 3.6; 
V_g = V_a; % No wind
g = 9.81;
K = 50000;
Ts = 0.01;
d = deg2rad(1.5); 

% System matrices
A = [-0.322,    0.052,   0.028, -1.12,   0.002; 
      0,        0,       1,     -0.001,  0; 
     -10.6,     0,      -2.87,   0.46,  -0.65; 
      6.87,     0,      -0.04,  -0.32,  -0.02; 
      0,        0,       0,      0,     -7.5];
B = [0, 0, 0, 0, 7.5]';
C = [1, 0, 0, 0, 0; 
     0, 1, 0, 0, 0; 
     0, 0, 1, 0, 0;
     0, 0, 0, 1, 0];

% Task defined constants for saturation
delta_a_max = abs(deg2rad(30.0));
e_phi_max   = abs(deg2rad(15.0));

% Transferfunction constants for Phi
a_phi_1 =  2.87;
a_phi_2 = -0.65;

% Task defined constant for Phi
zetta_phi = 0.707;

% Omega_n for Phi
w_n_phi = sqrt(abs(a_phi_2) * delta_a_max / e_phi_max);

% Phi PID constants 
k_p_phi = sign(a_phi_2) * delta_a_max / e_phi_max;
k_d_phi = (2 * zetta_phi * w_n_phi - a_phi_1) / a_phi_2;
% From root-locus analysis
% Should be in range [-pi, 0]
k_i_phi = -0.01 * 0; 

% User defined constants for Chi
W_chi = 10; 
zetta_chi = 2; 

% Omega_n for Chi
w_n_chi = w_n_phi / W_chi; 

% Chi PI constants 
k_p_chi = 2 * zetta_chi * w_n_chi * V_g / g; 
k_i_chi = w_n_chi^2 * V_g / g;

% Initial values
chi_0       = deg2rad(20);
beta_0      = 0;
p_0         = 0;
phi_0       = 0;
r_0         = 0; 
delta_a_0   = 0; 

% Indices
beta    = 1;
phi     = 2; 
p       = 3; 
r       = 4; 
delta_a = 5; 

%% Allocation 
% Time
t           = Ts * (0:(K-1));
% States
chi         = zeros(1, K); 
x           = zeros(5, K);

% Errors
e_chi       = zeros(1, K);
e_chi_int   = zeros(1, K);
e_phi       = zeros(1, K);
e_phi_unsat = zeros(1, K);
e_phi_int   = zeros(1, K);

% References
delta_a_ref         = zeros(1, K);
delta_a_ref_unsat   = zeros(1, K);
chi_ref             = zeros(1, K);
phi_ref             = zeros(1, K);

%% Initialization 
chi(1)  = chi_0;
x(:, 1) = [beta_0, phi_0, p_0, r_0, delta_a_0]';

% steps from 30 to 20 to 10 to 0 degs
chi_ref(1:K/4) = deg2rad(20); 
chi_ref(K/4:K/2) = deg2rad(10); 
chi_ref(K/2:3*K/4) = deg2rad(5); 

%% Simulation

for k = 1:K
    % Error in Chi
    e_chi(k) = chi_ref(k) - chi(k);
    
    % Phi_c, set based on PI-controller with error in Chi
    phi_ref(k) = k_i_chi * e_chi_int(k) + k_p_chi * e_chi(k); 
    
    % Error in Phi
    e_phi_unsat(k) = phi_ref(k) - x(phi, k);
    e_phi(k) = min(e_phi_max, max(-e_phi_max, e_phi_unsat(k))); % abs(e_phi) <= e_phi_max
    
    % Delta_a^c, set based on PID-controller with error in Phi
    delta_a_ref_unsat(k) = k_i_phi * e_phi_int(k) + k_p_phi * e_phi(k) - k_d_phi * x(p, k);
    delta_a_ref(k) = min(delta_a_max, max(-delta_a_max, delta_a_ref_unsat(k))); % abs(delta_a) <= delta_a_max
    
    if k < K
        % integrate states
        % x(:, k + 1) = euler2(A * x(:, k) + B * delta_a_ref(k), x(:, k), Ts);
        x(p, k + 1) = euler2(a_phi_2 * delta_a_ref(k) - a_phi_1 * x(p, k), x(p, k), Ts);
        x(phi, k + 1) = euler2(x(p, k), x(phi, k), Ts);
        chi(k + 1) = euler2((g / V_g) * (x(phi, k) + d), chi(k), Ts);
        % chi(k + 1) = euler2((g / V_g) * tan(x(phi, k) + d) * cos(x(beta, k)), chi(k), Ts); 
        
        % integrate errors
        if delta_a_ref(k) == delta_a_ref_unsat(k)
            e_chi_int(k + 1) = euler2(e_chi(k), e_chi_int(k), Ts);
        else
            e_chi_int(k + 1) = e_chi_int(k);
        end
        if e_phi(k) == e_phi_unsat(k)
            e_phi_int(k + 1) = euler2(e_phi(k), e_phi_int(k), Ts);
        else
            e_phi_int(k + 1) = e_phi_int(k);
        end
%         if k_i_chi ~= 0
%             e_chi_int(k + 1) = euler2((1/k_i_chi) * (e_phi(k) - e_phi_unsat(k)), e_chi_int(k + 1), Ts);
%         end
%         if k_i_phi ~= 0
%             e_phi_int(k + 1) = euler2((1/k_i_phi) * (delta_a_ref(k) - delta_a_ref_unsat(k)), e_phi_int(k + 1), Ts); 
%         end
    end
end

%% Plotting
if to_plot
    % Course (Chi)
    fig1 = figure(1); clf;
    plot(t, rad2deg(chi), t, rad2deg(chi_ref));
    legend('Course (Chi)', 'Course reference');
    ylabel('Course [deg]'); 
    xlabel('Time [s]');
    grid on; 
    
    if to_print
        set(fig1, 'Units', 'Inches');
        pos1 = get(fig1, 'Position');
        set(fig1, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos1(3), pos1(4)]);
        print(fig1, '3e_chi_course', '-depsc', '-r0');
    end
        
    % Aileron (delta_a)
    fig2 = figure(2); clf;
    plot(t, rad2deg(delta_a_ref), t, rad2deg(delta_a_max) * ones(1, K), t, - rad2deg(delta_a_max) * ones(1, K)); 
    legend('Aileron', 'Aileron positive saturation', 'Aileron negative saturation'); 
    ylabel('Aileron [deg]');
    xlabel('Time [s]');
    ylim([-rad2deg(delta_a_max) * 1.1, rad2deg(delta_a_max) * 1.1]); % increase view by 10% 
    grid on; 

    if to_print
        set(fig2, 'Units', 'Inches');
        pos1 = get(fig2, 'Position');
        set(fig2, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos1(3), pos1(4)]);
        print(fig2, '3e_delta_a_aileron', '-depsc', '-r0');
    end
end
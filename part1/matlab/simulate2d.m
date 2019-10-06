%% Clear
clear; close all;
to_plot = false;
%% Constants
% Global constants 
V_a = 580 / 3.6; 
V_g = V_a; % No wind
g = 9.81;
K = 50000;
Ts = 0.01;
d = deg2rad(1.5); 

% Task defined constants for saturation
delta_a_max = deg2rad(30.0);
e_phi_max = deg2rad(15.0);

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
chi_0 = deg2rad(20);
p_0   = 0;
phi_0 = 0;

%% Allocation 
% Time
t           = Ts * (0:(K-1));
% States
chi         = zeros(1, K); 
p           = zeros(1, K); 
phi         = zeros(1, K); 

% Errors
e_chi       = zeros(1, K);
e_chi_int   = zeros(1, K);
e_phi       = zeros(1, K);
e_phi_int   = zeros(1, K);

% References
delta_a     = zeros(1, K);
chi_ref     = zeros(1, K);
phi_ref     = zeros(1, K);

%% Initialization 
chi(1) = chi_0;
p(1) = p_0;
phi(1) = phi_0;

% steps from 30 to 20 to 10 to 0 degs
chi_ref(1:K/4) = deg2rad(20); 
chi_ref(K/4:K/2) = deg2rad(10); 
chi_ref(K/2:3*K/4) = deg2rad(5); 

%% Simulation

for k = 1:K
    e_chi(k) = chi_ref(k) - chi(k);
    
    phi_ref(k) = k_i_chi * e_chi_int(k) + k_p_chi * e_chi(k); 
    
    e_phi(k) = phi_ref(k) - phi(k);
    e_phi(k) = min(e_phi_max, max(-e_phi_max, e_phi(k))); % abs(e_phi) <= e_phi_max
    
    delta_a(k) = k_i_phi * e_phi_int(k) + k_p_phi * e_phi(k) - k_d_phi * p(k);
    delta_a(k) = min(delta_a_max, max(-delta_a_max, delta_a(k))); % abs(delta_a) <= delta_a_max
    
    if k < K
        % integrate states
        % a_phi_2 / (s + a_phi_1) -> a_phi_2 * exp(-a_phi_1 * t)
        % a_phi_2 * exp(-a_phi_2 * t(k)) * delta_a(k);
        % euler2(a_phi_2 * delta_a(k) - a_phi_1 * p(k), p(k), Ts);
        p(k + 1) = euler2(a_phi_2 * delta_a(k) - a_phi_1 * p(k), p(k), Ts);
        phi(k + 1) = euler2(p(k), phi(k), Ts);
        chi(k + 1) = euler2((g / V_g) * (phi(k) + d), chi(k), Ts);
%         p(k + 1) = p(k) + (delta_a(k) * a_phi_2 - p(k) * a_phi_1) * Ts;
%         phi(k + 1) = phi(k) + p(k) * Ts; 
%         chi(k + 1) = chi(k) + (g/V_g) * (phi(k) + d) * Ts;
%         
%         e_chi_int(k + 1) = e_chi_int(k) + e_chi(k) * Ts;
%         e_phi_int(k + 1) = e_phi_int(k) + e_phi(k) * Ts;
        
        % integrate errors
        e_chi_int(k + 1) = euler2(e_chi(k), e_chi_int(k), Ts);
        e_phi_int(k + 1) = euler2(e_phi(k), e_phi_int(k), Ts);
    end
end

%% Plotting
% Course (Chi)
if to_plot
    fig1 = figure(1); clf;
    plot(t, rad2deg(chi), t, rad2deg(chi_ref));
    legend('Course (Chi)', 'Course reference');
    ylabel('Course [deg]'); 
    xlabel('Time [s]');
    grid on; 

    set(fig1, 'Units', 'Inches');
    pos1 = get(fig1, 'Position');
    set(fig1, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos1(3), pos1(4)]);
    print(fig1, '2d_chi_course', '-depsc', '-r0');

    % Aileron (delta_a)
    fig2 = figure(2); clf;
    plot(t, rad2deg(delta_a), t, rad2deg(delta_a_max) * ones(1, K), t, - rad2deg(delta_a_max) * ones(1, K)); 
    legend('Aileron', 'Aileron positive saturation', 'Aileron negative saturation'); 
    ylabel('Aileron [deg]');
    xlabel('Time [s]');
    ylim([-rad2deg(delta_a_max) * 1.1, rad2deg(delta_a_max) * 1.1]); % increase view by 10% 
    grid on; 

    set(fig2, 'Units', 'Inches');
    pos1 = get(fig2, 'Position');
    set(fig2, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos1(3), pos1(4)]);
    print(fig2, '2d_delta_a_aileron', '-depsc', '-r0');
end

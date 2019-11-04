%% Clear
clear; close all;
to_plot = true;
to_print = true;

only_kalman = true;
sensor_failiure = true;

if only_kalman
    plot_folder = '3f';
else
    plot_folder = '3e';
end

if sensor_failiure
    plot_folder = '3g';
    if ~only_kalman
        to_print = false;
    end
end

base_folder = 'part2\figures\';
%% Constants
% Global constants 
rng('default');
V_a = 580 / 3.6; 
V_g = V_a; % No wind
g = 9.81;
K = 500000;
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

% Measurement
C_m = [0, 0, 1, 0, 0; 
       0, 0, 0, 1, 0]; 
 
% % Kalman matrices
% A_k = [-0.322,  0.052,  0.028,  -1.12; 
%         0,      0,      1,      -0.001; 
%        -10.6,   0,     -2.87,    0.46; 
%         6.87,   0,     -0.04,   -0.32]; 
% B_k = [0.002, 0, -0.65, -0.02]'; 
% C_k = [0, 0, 1, 0; 
%        0, 0, 0, 1]; 
% E_k = eye(4); 

% Noise
% Process noise
Q = Ts * 1e-6 * diag([0.001, 1, 100, 10, 0]);

% Measurement noise
R = deg2rad(diag([0.2, 0.2])).^2;

% Kalman noise
Q_k = 1.2 * Q(:, :);
R_k = R(:, :);
 
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

% Simulator initial values
chi_0       = deg2rad(20);
beta_0      = 0;
p_0         = 0;
phi_0       = 0;
r_0         = 0; 
delta_a_0   = 0; 

% Estimator initial values
chi_bar_0   = deg2rad(15); 
beta_bar_0  = 0;
p_bar_0     = 0;
phi_bar_0   = 0;
r_bar_0     = 0; 
P_bar_0     = 1e-5 * eye(4); % zeros(4, 4); 

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
% Estimator
x_hat       = zeros(4, K);
P_hat       = zeros(4, 4, K); 
x_bar       = zeros(4, K);
P_bar       = zeros(4, 4, K); 
% Measurement
z           = zeros(2, K); 

% Errors
e_chi               = zeros(1, K);
e_chi_int           = zeros(1, K);
e_phi               = zeros(1, K);
e_phi_unsat         = zeros(1, K);
e_phi_saturated     = false(1, K); 
e_phi_int           = zeros(1, K);

% References
delta_a_ref         = zeros(1, K);
delta_a_ref_unsat   = zeros(1, K);
delta_a_saturated   = false(1, K); 
chi_ref             = zeros(1, K);
phi_ref             = zeros(1, K);

%% Initialization 
chi(1)          = chi_0;
x(:, 1)         = [beta_0, phi_0, p_0, r_0, delta_a_0]';
x_bar(:, 1)     = [beta_bar_0, phi_bar_0, p_bar_0, r_bar_0]'; 
P_bar(:, :, 1)  = P_bar_0; 

% steps from 30 to 20 to 10 to 0 degs
chi_ref(1:K/4)      = deg2rad(20); 
chi_ref(K/4:K/2)    = deg2rad(10); 
chi_ref(K/2:3*K/4)  = deg2rad(5); 

%% Simulation
for k = 1:K
    % Sensor failiure
    if sensor_failiure
        if k > K/2
            R_k = 1e100 * R_k; 
            sensor_failiure = false;
        end
    end
    
    % Measurement
    z(:, k) = C_m * x(:, k) + mvnrnd(zeros(1, 2), R.^(0.5))'; 
    
    % Error in Chi
    e_chi(k) = chi_ref(k) - chi(k); % TODO: Find out: Should chi still be used like this? 
    
    % Phi_c, set based on PI-controller with error in Chi
    phi_ref(k) = k_i_chi * e_chi_int(k) + k_p_chi * e_chi(k); 
    
    % Error in Phi
    % e_phi_unsat(k) = phi_ref(k) - x_hat(phi, k);
    e_phi_unsat(k) = phi_ref(k) - x_bar(phi, k);
    [e_phi(k), e_phi_saturated(k)] = saturate(e_phi_unsat(k), e_phi_max); 
    % e_phi(k) = min(e_phi_max, max(-e_phi_max, e_phi_unsat(k))); % abs(e_phi) <= e_phi_max
    
    % Delta_a^c, set based on PID-controller with error in Phi
    % delta_a_ref_unsat(k) = k_i_phi * e_phi_int(k) + k_p_phi * e_phi(k) - k_d_phi * x_hat(p, k);
    if only_kalman
        delta_a_ref_unsat(k) = k_i_phi * e_phi_int(k) + k_p_phi * e_phi(k) - k_d_phi * x_bar(p, k);
    else
        delta_a_ref_unsat(k) = k_i_phi * e_phi_int(k) + k_p_phi * e_phi(k) - k_d_phi * z(1, k);
    end
    [delta_a_ref(k), delta_a_saturated(k)] = saturate(delta_a_ref_unsat(k), delta_a_max); 
    % delta_a_ref(k) = min(delta_a_max, max(-delta_a_max, delta_a_ref_unsat(k))); % abs(delta_a) <= delta_a_max
    
    if k < K
        % Simulate actual system based on calculated input
        x(:, k + 1) = euler2(A * x(:, k) + B * delta_a_ref(k) + mvnrnd(zeros(1, 5), Q.^(0.5))', x(:, k), Ts);
        chi(k + 1) = euler2((g / V_g) * tan(x_bar(phi, k) + d) * cos(x_bar(beta, k)), chi(k), Ts); 
        
        
        % Kalman filter
        [x_bar_next, P_bar_next, x_hat_, P_hat_] = ...
            KalmanFilter(x_bar(:, k), P_bar(:, :, k), delta_a_ref(k), z(:, k), Ts, R_k, Q_k(beta:r, beta:r));
        x_hat(:, k) = x_hat_; 
        P_hat(:, :, k) = P_hat_; 
        x_bar(:, k + 1) = x_bar_next; 
        P_bar(:, :, k + 1) = P_bar_next; 
        % [x_post, P_post, x_pred, P_pred, v_inno, S_inno] = ...
            % KalmanFilter(x_hat(beta:r, k), P_hat(beta:r, beta:r, k), delta_a_ref(k), z(:, k), A_k, B_k, C_k, Q, R);
        % x_hat(:, k + 1) = x_post; 
        % P_hat(:, :, k + 1) = P_post; 
        
        % Integrate errors
        if e_phi_saturated(k)
            e_chi_int(k + 1) = e_chi_int(k);
        else
            e_chi_int(k + 1) = euler2(e_chi(k), e_chi_int(k), Ts);
        end
        if delta_a_saturated(k)
            e_phi_int(k + 1) = e_phi_int(k);
        else
            e_phi_int(k + 1) = euler2(e_phi(k), e_phi_int(k), Ts);
        end
    end
end

%% Plotting
% plot_folder = '3e';
if to_plot
    base_folder
    plot_folder
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
        print(fig1, strcat(base_folder, plot_folder, '/', 'chi_course'), '-depsc', '-r0');
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
        print(fig2, strcat(base_folder, plot_folder, '/delta_a_aileron'), '-depsc', '-r0');
    end
    
    % Estimated and true sideslip angle
    fig3 = figure(3); clf; 
    plot(t, rad2deg(x_bar(beta, :)), t, rad2deg(x(beta, :))); 
    legend('Estimated sideslip', 'True sideslip'); 
    ylabel('Sideslip angle [deg]');
    xlabel('Time [s]');
    grid on; 
    
    if to_print
        set(fig3, 'Units', 'Inches'); 
        pos1 = get(fig3, 'Position'); 
        set(fig3, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos1(3), pos1(4)]);
        print(fig3, strcat(base_folder, plot_folder, '/', 'beta_sideslip'), '-depsc', '-r0');
    end
    
    % Estimated and true roll angle
    fig4 = figure(4); clf; 
    if sensor_failiure
        plot(t, rad2deg(x_bar(phi, :)), t, rad2deg(x(phi, :)), t, abs(rad2deg(x_bar(phi, :)) - rad2deg(x(phi, :)))); 
        legend('Estimated roll', 'True roll', 'Difference'); 
    else
        plot(t, rad2deg(x_bar(phi, :)), t, rad2deg(x(phi, :))); 
        legend('Estimated roll', 'True roll'); 
    end
    % plot(t, rad2deg(x_bar(phi, :)), t, rad2deg(x(phi, :))); 
    % legend('Estimated roll', 'True roll'); 
    ylabel('Roll angle [deg]');
    xlabel('Time [s]');
    grid on; 
    
    if to_print
        set(fig4, 'Units', 'Inches'); 
        pos1 = get(fig4, 'Position'); 
        set(fig4, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos1(3), pos1(4)]);
        print(fig4, strcat(base_folder, plot_folder, '/', 'roll_phi'), '-depsc', '-r0');
    end
    
    % Estimated, noisy and true roll rate
    fig5 = figure(5); clf; 
    plot(t, rad2deg(z(1, :)), t, rad2deg(x_bar(p, :)), t, rad2deg(x(p, :))); 
    legend('Noisy roll rate', 'Estimated roll rate', 'True roll rate'); 
    ylabel('Roll rate angle [deg/s]');
    xlabel('Time [s]');
    grid on; 
    
    if to_print
        set(fig5, 'Units', 'Inches'); 
        pos1 = get(fig5, 'Position'); 
        set(fig5, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos1(3), pos1(4)]);
        print(fig5, strcat(base_folder, plot_folder, '/', 'roll_rate_p'), '-depsc', '-r0');
    end
    
    % Estimated, noisy and true yaw rate
    fig6 = figure(6); clf; 
    plot(t, rad2deg(z(2, :)), t, rad2deg(x_bar(r, :)), t, rad2deg(x(r, :))); 
    legend('Noisy yaw rate', 'Estimated yaw rate', 'True yaw rate'); 
    ylabel('Yaw rate angle [deg/s]');
    xlabel('Time [s]');
    grid on; 
    
    if to_print
        set(fig6, 'Units', 'Inches'); 
        pos1 = get(fig6, 'Position'); 
        set(fig6, 'PaperPositionMode', 'Auto', 'PaperUnits', 'Inches', 'PaperSize', [pos1(3), pos1(4)]);
        print(fig6, strcat(base_folder, plot_folder, '/', 'yaw_rate_r'), '-depsc', '-r0');
    end
end



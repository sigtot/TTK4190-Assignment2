function [x_bar_next, P_bar_next, x_hat, P_hat] = ...
    KalmanFilter(x_bar_prev, P_bar_prev, u_prev, z, Ts, R, Q)
%KalmanFilter Summary of this function goes here
%   Detailed explanation goes here
    A_k = [-0.322,  0.052,  0.028,  -1.12; 
            0,      0,      1,      -0.001; 
           -10.6,   0,     -2.87,    0.46; 
            6.87,   0,     -0.04,   -0.32]; 
    B_k = [0.002, 0, -0.65, -0.02]'; 
    C_k = [0, 0, 1, 0; 
           0, 0, 0, 1]; 
    E_k = eye(4);
    
    A_d = @(Ts) eye(size(A_k)) + Ts * A_k;
    B_d = @(Ts) Ts * B_k;
    C_d = C_k; 
    E_d = E_k; % multiply with Ts?

    PHI   = A_d(Ts); 
    DELTA = B_d(Ts); 
    GAMMA = E_d; 
    
    H = C_d; 
    
    K = P_bar_prev * H' / (H * P_bar_prev * H' + R);
    x_hat = x_bar_prev + K * (z - H * x_bar_prev); 
    P_hat = (eye(4) - K * H) * P_bar_prev * (eye(4) - K * H)' ...
        + K * R * K'; 
    
    x_bar_next = PHI * x_hat + DELTA * u_prev;
    P_bar_next = PHI * P_hat * PHI' + GAMMA * Q * GAMMA'; 
end
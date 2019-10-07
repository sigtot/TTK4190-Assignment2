function [x_post, P_post, x_pred, P_pred, v_inno, S_inno] = ...
    KalmanFilter(x_prev, P_prev, u_prev, z, A, B, C, Q, R)
%KalmanFilter Summary of this function goes here
%   Detailed explanation goes here
    x_pred = A * x_prev + B * u_prev; 
    P_pred = A * P_prev * A' + Q; 
    
    z_pred = C * x_pred;
    
    v_inno = z - z_pred;
    S_inno = C * P_pred * C' + R;
    
    W_kalm = P_pred * C' / S_inno; 
    
    x_post = x_pred + W_kalm * v_inno; 
    c = W_kalm * C; 
    P_post = (eye(size(c)) - c) * P_pred;
end


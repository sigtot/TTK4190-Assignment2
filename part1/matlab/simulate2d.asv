%% Constants
V_a = 580; 
V_g = 580 / 3.6;
g = 9.81;
K = 100;
Ts = 0.1;
d = deg2rad(1.5); 

delta_a_max = deg2rad(30.0);
e_phi_max = deg2rad(15.0);

a_phi_1 =  2.87;
a_phi_2 = -0.65;

zetta_phi = 0.707;

w_n_phi = sqrt(abs(a_phi_2) * delta_a_max / e_phi_max);

k_p_phi = sign(a_phi_2) * delta_a_max / e_phi_max;
k_d_phi = (2 * zetta_phi * w_n_phi - a_phi_1) / a_phi_2;
k_i_phi = 0; % from root-locus analysis

W_chi = 10; 
zetta_chi = 2; 

w_n_chi = w_n_phi / W_chi; 

k_p_chi = 2 * zetta_chi * w_n_chi * V_g / g; 
k_i_chi = w_n_chi^2 * V_g / g;

%% Allocation 
chi         = zeros(1, K); 
p           = zeros(1, K); 
phi         = zeros(1, K); 
e_chi_int   = zeros(1, K);
e_phi_int   = zeros(1, K);

delta_a     = zeros(1, K);
chi_ref     = zeros(1, K);

%% Initialization 
chi(1) = chi_0;
p(1) = p_0;
phi(1) = phi_0;

%% Simulation


for k = 1:K
    
end
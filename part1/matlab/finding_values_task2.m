a_phi_1 =  2.87;
a_phi_2 = -0.65;

delta_a_max = deg2rad(30.0);
e_phi_max = deg2rad(15.0);

zetta_phi = 0.707;

w_n_phi = sqrt(abs(a_phi_2) * delta_a_max / e_phi_max);

k_p_phi = sign(a_phi_2) * delta_a_max / e_phi_max;
k_d_phi = (2 * zetta_phi * w_n_phi - a_phi_1) / a_phi_2;

sys = tf([a_phi_2], [1, a_phi_1 + a_phi_2 * k_d_phi, a_phi_2 * k_p_phi, 0]);
k_i_phi = -100:100;

rlocus(sys, k_i_phi);


%% solver
% time
t_i = 0;
s_max = 0.005;
t_f = 0;
% gravitation
g = [0,-9.80665,0];

%% scapulohumeral rhythm
% translation
x_ej = -0.05;
y_ej = 0.06;
z_ej = -0.03;
% euler angles
psi_ej = 0; % rot y
mu_ej = 0; % rot x
phi_ej = 0; % rot y
% ellipsoidal joint
r_ej = 0.087;
h_ej = 0.15;

k_ej = 100;
b_ej = 10;

h_min = -0.9*h_ej;
h_max = +0.9*h_ej;
% time delay
T_d = 1e-20;
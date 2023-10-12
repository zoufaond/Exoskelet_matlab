%% solver
% time
t_i = 0;
s_max = 0.005;
t_f = 0;
% gravitation
g = [0,-9.80665,0];
f_koef = 1e-3;

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
r_ej = 0.0821;
h_ej = 0.15;

k_ej = 100;
b_ej = 10;

h_min = -0.9*h_ej;
h_max = +0.9*h_ej;
% time delay
T_d = 1e-20;


%%% Muscles init %%%
Glenohumeral = readmatrix('Glenohumeral_muscles.xlsx');
Glenohumeral_F0M = Glenohumeral(2:end,2);
Glenohumeral_l0  = Glenohumeral(2:end,3)*4;

Scapulothoracic = readmatrix('Scapulothoracic_muscles.xlsx');
Scapulothoracic_F0M = Scapulothoracic(2:end,2);
Scapulothoracic_l0  = Scapulothoracic(2:end,3)*4;
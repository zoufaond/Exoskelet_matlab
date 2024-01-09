clc; clear all; close all;

%% solver
% time
t_i = 0;
s_max = 0.005;
t_f = 0.5;
% gravitation
g = [0;0;0]; % [0,-9.80665,0];

%% scapulohumeral rhythm
% translation
x_ej = -0.03; %zmena
y_ej = -0.1; % %zmena
z_ej = 0.05; %zmena
% euler angles
psi_ej = 0; % rot y
mu_ej = 0; % rot x
phi_ej = 0; % rot y
% ellipsoidal joint
r_ej = 0.088; % 0.087;
h_ej = 0.19; % 0.15;
h_min = -0.95*h_ej;
h_max = +0.95*h_ej;
% time delay
T_d = 1e-20;
% binding springs
k_S = 500;
b_S = 50;

%%% Muscles init %%%
Glenohumeral = readmatrix('Glenohumeral_muscles.xlsx');
Glenohumeral_F0M = Glenohumeral(2:end,2)*0.01;
Glenohumeral_l0  = Glenohumeral(2:end,3);

Scapulothoracic = readmatrix('Scapulothoracic_muscles.xlsx');
Scapulothoracic_F0M = Scapulothoracic(2:end,2)*0.01;
Scapulothoracic_l0  = Scapulothoracic(2:end,3);

GH = readtable('Glenohumeral_muscles.xlsx');
Glenohumeral_F0M = str2double(table2array(GH(:,2)))';
Glenohumeral_l0 = str2double(table2array(GH(:,3)))';
ST = readtable('Scapulothoracic_muscles.xlsx');
Scapulothoracic_F0M = str2double(table2array(ST(:,2)))';
Scapulothoracic_l0 = str2double(table2array(ST(:,3)))';

%real konstanty
konst = 1;
cI = [5.75482e-05, 6.05521e-05, 1.05345e-05,-1.27405e-05, 1.65861e-05, 4.50136e-06];
ccom = [-0.0110972, 0.00637508, 0.0541825];
sI = [0.00020065, 0.000202939, 0.000227237,-4.55661e-05, -7.53511e-05, -7.42694e-05];
scom = [-0.0512948, -0.0367122, -0.0430756];
hI = [0.00326632, 4.34807e-05, 0.00325004,-0.000136456, -1.30958e-06, -1.14964e-05];
hcom = [0.000389382, -0.116338, -0.00231072];
m = [0.0370005,0.12759,0.26363];
T_t = [0.006325 0.00693 0.025465];
T_c = [-0.01433 0.02007 0.135535];
T_s = [-0.00955 -0.034 0.009];
c = 0;
k = 0;
g = 9.8066;

p = [c,k,g];

%% scapulohumeral rhythm
%kontaktni body na scapule
cont1 = [-0.083529 -0.01368 -0.098954]-scom;
cont2 = [-0.105456 -0.125076 -0.075251]-scom;
% translation
x_ej = -0.03; %zmena
y_ej = -0.1; % %zmena
z_ej = 0.05; %zmena
mx = x_ej;
my = y_ej;
mz = z_ej;
m_el = [x_ej,y_ej,z_ej];
% euler angles
psi_ej = 0; % rot y
mu_ej = 0; % rot x
phi_ej = 0; % rot y
% ellipsoidal joint
r_ej = 0.088; % 0.087;
h_ej = 0.15; % 0.15;
ax = r_ej;
ay = h_ej;
az = r_ej;
a_el = [r_ej,h_ej,r_ej];
h_min = -0.95*h_ej;
h_max = +0.95*h_ej;
% time delay
T_d = 1e-20;
% binding springs
k_S = 500;
b_S = 50;
% contact force Chadwick
k_contact = 0;
eps = 1e-3;
cont_params = [k_contact, eps];

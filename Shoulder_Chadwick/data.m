g = 9.81;
SC_P_offset = [0.0014 -0.0152 0.0028];
AC_P_offset = [0.2077 -0.0001 0.0282];
GH_P_offset = [0.002 -0.023 -0.0439];

SC_yzx = [-33.4 5.002 32.914]*pi/180;
AC_yzx = [45.571 0.458 -12.062]*pi/180;
GH_yzy = [43.864 -0.218 -34.659]*pi/180;

clavicle_com = [0.0983 0.0167 0.0042];
scapula_com = [-0.053 -0.0229 -0.0247];
humerus_com = [0.0064 -0.0776 -0.0036];
konst = 1e1;
IU = [6.4e-06 2.63e-05 2.43e-05 -6.7e-06 -9.5e-06 2.9e-06]*konst;
IM = [0.001 0.001 0.001 0 0 0]*konst;
IL = [0.0132 0.001988 0.0132 0 0 0]*konst;
m = [0.156,0.7054,2.0519];

rigid0P = SC_P_offset;
rigid1C = -clavicle_com;
rigid1P = -clavicle_com+AC_P_offset;
rigid2C = -scapula_com;
rigid2P = -scapula_com+GH_P_offset;
rigid3C = -humerus_com;

c = 1;
k = 0;

%% contact_data
%kontaktni body na scapule
cont1 = [-0.124743  0.003165 -0.007732]-scapula_com;
cont2 = [-0.109446 -0.130736 -0.005218]-scapula_com;
% translation
x_ej = -0.0; %zmena
y_ej = -0.1486; % %zmena
z_ej = 0.0591; %zmena
mx = x_ej;
my = y_ej;
mz = z_ej;
m_el = [x_ej,y_ej,z_ej];
ax = 0.147;
ay = 0.2079;
az = 0.0944;
% contact force Chadwick
k_contact = 100000;
eps = 1e-4;
cont_params = [mx,my,mz,ax,ay,az,cont1(1),cont1(2),cont1(3),cont2(1),cont2(2),cont2(3),k_contact, eps];
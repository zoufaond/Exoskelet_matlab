g = 9.81;
SC_P_offset = [0.0014 -0.0152 0.0028];
AC_P_offset = [0.2077 -0.0001 0.0282];
GH_P_offset = [0.002 -0.023 -0.0439];
EL_P_offset = [0.0058 -0.2907 0.0049];
EL_P_offset_rot = [0 0 0.32318];
EL_rot_axis = [0.969 -0.247 0];
[ELx_ax,ELx_angle] = vecs2ax_angle(EL_rot_axis,[0,0,1]);
PS_P_offset = [0.0067234 -0.00931306 -0.0068107];
PSY_rot_axis = [0.182 0.98227 -0.044946];
[PSy_ax,PSy_angle] = vecs2ax_angle(PSY_rot_axis,[0,0,1]);
hand_P_offset = [-0.0333 -0.2342 0.0092];


SC_yzx = [-21.784 6.303 0]*pi/180;
AC_yzx = [46.295 4.899 -1.180]*pi/180;
GH_yzy = [0 40 0]*pi/180;
EL_x = 150*pi/180;
PS_y = 100*pi/180;

clavicle_com = [0.0983 0.0167 0.0042];
scapula_com = [-0.053 -0.0229 -0.0247];
humerus_com = [0.0064 -0.0776 -0.0036];
ulna_com = [-0.0199622 -0.0594577 0.0101954];
radius_com = [-0.00736417 -0.118863 0.00170825];
hand_com = [0.0006 -0.0905 -0.0365];

konst = 1e1;
IU = [6.4e-06 2.63e-05 2.43e-05 0 0 0]*konst;
IM = [0.001 0.001 0.001 0 0 0]*konst;
IL = [0.0132 0.001988 0.0132 0 0 0]*konst;
m = [2.0519,0.5464,0.5464,0.525];

rigid0P = SC_P_offset;
rigid1C = -clavicle_com;
rigid1P = -clavicle_com+AC_P_offset;
rigid2C = -scapula_com;
rigid2P = -scapula_com+GH_P_offset;
rigid3C = -humerus_com;
rigid3P = -humerus_com+EL_P_offset;
rigid4C = -ulna_com;
rigid4P = -ulna_com+PS_P_offset;
rigid5C = -radius_com;
rigid5P = -radius_com+hand_P_offset;
rigid6C = -hand_com;

c = 10;
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

function [axis,angle] = vecs2ax_angle(rot_axis,base_ax)
    perp_axis = cross(base_ax,rot_axis);
    axis = perp_axis/norm(perp_axis);
    angle = acos(dot(rot_axis,base_ax));
end
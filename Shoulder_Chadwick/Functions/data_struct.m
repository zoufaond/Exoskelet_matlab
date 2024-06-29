addpath('EoMs_derivation\')
addpath('Functions\')
load('das3_simplified.mat')
load('data_model.mat')
muscles = model_simpl.muscles;
for i=1:length(muscles)
    fmax_vec(i) = muscles{i}.fmax*4;
    lceopt_vec(i) = muscles{i}.lceopt+0;
    lslack_vec(i) = muscles{i}.lslack;
end

%initial conditions
SC_yzx = [-21.784 6.303 0]*pi/180;
AC_yzx = [46.295 4.899 -1.180]*pi/180;
GH_yzy = [0 5 0]*pi/180;
EL_x = 0*pi/180;
PS_y = 5*pi/180;
initCond = [SC_yzx,AC_yzx,GH_yzy,EL_x,zeros(1,10)]';
activ = zeros(35,1);

%global constants
model.g = 9.81;
model.c = 0.5;
model.k = 0;

%offset translations and rotations in each joints (measured from parents)
model.SC_P_offset = [0.0014 -0.0152 0.0028];
model.AC_P_offset = [0.2077 -0.0001 0.0282];
model.GH_P_offset = [0.002 -0.023 -0.0439];
model.EL_P_offset = [0.0058 -0.2907 0.0049];
model.EL_P_offset_rot = [0 0 0.32318];
model.EL_rot_axis = [0.969 -0.247 0];
[model.ELx_ax,model.ELx_angle] = vecs2ax_angle(model.EL_rot_axis,[0,0,1]);
model.PS_P_offset = [0.0067234 -0.00931306 -0.0068107];
model.PSY_rot_axis = [0.182 0.98227 -0.044946];
[model.PSy_ax,model.PSy_angle] = vecs2ax_angle(model.PSY_rot_axis,[0,0,1]);
model.hand_P_offset = [-0.0333 -0.2342 0.0092];

% center of mass of segments
model.clavicula_com = [0.0983 0.0167 0.0042];
model.scapula_com = [-0.053 -0.0229 -0.0247];
model.humerus_com = [0.0064 -0.0776 -0.0036];
model.ulna_com = [-0.0199622 -0.0594577 0.0101954];
model.radius_com = [-0.00736417 -0.118863 0.00170825];
model.hand_com = [0.0006 -0.0905 -0.0365];

% segments inertias
model.I_clavicula = [6.4e-06 2.63e-05 2.43e-05 0 0 0];
model.I_scapula = [0.001 0.001 0.001 0 0 0];
model.I_humerus = [0.0132 0.001988 0.0132 0 0 0];
model.I_ulna = [0.0030585 0.00045325 0.0030585 0 0 0];
model.I_radius = [0.0030585 0.00045325 0.0030585 0 0 0];
model.I_hand = [0.0006387 0.0001904 0.0006387 0 0 0];

% segment masses
model.clavicula_mass = 0.156;
model.scapula_mass = 0.7054;
model.humerus_mass = 2.0519;
model.ulna_mass = 0.5464;
model.radius_mass = 0.5464;
model.hand_mass = 0.525;

% recalculated transformation w.r.t. segments coms
model.rigid0P = model.SC_P_offset;
model.rigid1C = -model.clavicula_com;
model.rigid1P = -model.clavicula_com+model.AC_P_offset;
model.rigid2C = -model.scapula_com;
model.rigid2P = -model.scapula_com+model.GH_P_offset;
model.rigid3C = -model.humerus_com;
model.rigid3P = -model.humerus_com+model.EL_P_offset;
model.rigid4C = -model.ulna_com;
model.rigid4P = -model.ulna_com+model.PS_P_offset;
model.rigid5C = -model.radius_com;
model.rigid5P = -model.radius_com+model.hand_P_offset;
model.rigid6C = -model.hand_com;

%% contact_data
%kontaktni body na scapule
model.contTS = [-0.124743  0.003165 -0.007732]-model.scapula_com;
model.contAI = [-0.109446 -0.130736 -0.005218]-model.scapula_com;
% translation
model.elips_trans = [0 -0.1486 0.0591];
model.elips_dim = [0.147 0.2079 0.0944];
% contact force Chadwick
model.k_contact_in = 20e5;
model.eps_in = 1e-5;
model.k_contact_out = 20;
model.eps_out = 1e-1;

%%
save('data_model.mat','model')
%%

function [axis,angle] = vecs2ax_angle(rot_axis,base_ax)
    perp_axis = cross(base_ax,rot_axis);
    axis = perp_axis/norm(perp_axis);
    angle = acos(dot(rot_axis,base_ax));
end
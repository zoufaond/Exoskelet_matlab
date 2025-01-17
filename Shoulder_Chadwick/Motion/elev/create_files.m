addpath ..\..\Functions\
addpath ..\
clearvars
motion_file = readmatrix('scabduction.txt');
% data_model = load('mot_struct_rigged.mat');
% 
% coords_struct_clavmod.tout = data_model.coords_struct.tout;
% coords_struct_clavmod.data = data_model.coords_struct.mot_euler_mod;
% coords_struct_clavmod.data(:,8) = coords_struct_clavmod.data(:,8) - coords_struct_clavmod.data(2,8);
% data2mot(coords_struct_clavmod,'abd_rigged.mot','euler', 'struct');


% initPosOpt = data_model.params.InitPosOptEul;
mot_unpoly.tout = motion_file(:,1)-motion_file(1,1);
mot_unpoly.data = motion_file(:,8:end-1)*pi/180;
% 
res = motion_polyfit(motion_file,'no plot_polyfits');
res.data(:,8) = res.data(:,8) - 0.3 * (0.94 - res.data(:,8));
data2mot(res,'scabduction_GL.mot','euler', 'struct');
% data2mot(mot_unpoly,'steering_not_polyfitted.mot','euler', 'struct')
% 
% coords_struct.tout = motion_file(:,1)-motion_file(1,1);
% coords_struct.mot_euler = res.data;
% % rigging = coords_struct.mot_euler(15,8);
% % coords_struct.mot_euler(:,8) = coords_struct.mot_euler(:,8);
% 
% 
% % coords_struct.mot_euler_IC(:,8) = coords_struct.mot_euler_IC(:,8);
% 
% coords_struct.mot_quaternion = eul2quat_motion(coords_struct.mot_euler);
% 
% coords_struct_clavmod.tout = motion_file(:,1)-motion_file(1,1);
% coords_struct_clavmod.data = change_clavx(coords_struct.mot_euler);
% data2mot(coords_struct_clavmod,'steering_clavmod.mot','euler', 'struct');
% 
% coords_struct.mot_euler_IC = coords_struct_clavmod.data(1,:);
% coords_struct.mot_quaternion_IC = eul2quat_motion(coords_struct.mot_euler_IC);
% 
% coords_struct.mot_euler_mod = coords_struct_clavmod.data;
% coords_struct.mot_quaternion_mod = eul2quat_motion(coords_struct_clavmod.data);
% 
% save('steering_struct.mat','coords_struct')

addpath ..\Functions\
addpath ..\
clearvars
motion_file = readmatrix('steering.txt');
data_model = load('data_model.mat');
initPosOpt = data_model.params.InitPosOptEul;
mot_unpoly.tout = motion_file(:,1)-motion_file(1,1);
mot_unpoly.data = motion_file(:,8:end-1)*pi/180;

res = motion_polyfit(motion_file,'');
data2mot(res,'steering_polyfitted.mot','euler', 'struct');
data2mot(mot_unpoly,'steering_not_polyfitted.mot','euler', 'struct')

coords_struct.tout = motion_file(:,1)-motion_file(1,1);
coords_struct.mot_euler = res.data;
% rigging = coords_struct.mot_euler(15,8);
% coords_struct.mot_euler(:,8) = coords_struct.mot_euler(:,8);


% coords_struct.mot_euler_IC(:,8) = coords_struct.mot_euler_IC(:,8);

coords_struct.mot_quaternion = eul2quat_motion(coords_struct.mot_euler);

coords_struct_clavmod.tout = motion_file(:,1)-motion_file(1,1);
coords_struct_clavmod.data = change_clavx(coords_struct.mot_euler);
data2mot(coords_struct_clavmod,'steering_clavmod.mot','euler', 'struct');

coords_struct.mot_euler_IC = coords_struct_clavmod.data(1,:);
coords_struct.mot_quaternion_IC = eul2quat_motion(coords_struct.mot_euler_IC);

coords_struct.mot_euler_mod = coords_struct_clavmod.data;
coords_struct.mot_quaternion_mod = eul2quat_motion(coords_struct_clavmod.data);

save('steering_struct.mat','coords_struct')

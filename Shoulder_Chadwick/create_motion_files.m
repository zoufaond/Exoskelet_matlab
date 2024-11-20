clearvars
motion = 'Motions\Abduction\abduction.mot';
[motion_path,motion_name,extension] = fileparts(motion);
osim_file = 'das3.osim';
mydir = 'Polyfiles';
musclepoly_file = 'musclepoly';

das3_polynomials(osim_file,mydir,motion,musclepoly_file);
model = das3_readosim(osim_file,[motion_path,'\',mydir '\' musclepoly_file]);
save([motion_path,'\','OS_model'],'model')
motion2struct(motion)



function motion2struct(motion)

[motion_path,motion_name,extension] = fileparts(motion);

data_motion = readmatrix([motion_path,'\',motion_name,'.txt']);
time = data_motion(:,1);
angles = data_motion(:,5:end-1)*pi/180;
quats = eul2quat_motion(angles);

mot_struct.euler = angles;
mot_struct.quat = quats;
mot_struct.time = time;

save([motion_path,'\',motion_name],'mot_struct');
end

function res = eul2quat_motion(angles)

SC_EUL = angles(:,1:3);
AC_EUL = angles(:,4:6);
GH_EUL = angles(:,7:9);

SC_Q = eul2quat(SC_EUL,'YZX');
AC_Q = eul2quat(AC_EUL,'YZX');
GH_Q = eul2quat(GH_EUL,'YZY');

EL_x = angles(:,10);

res = [SC_Q,AC_Q,GH_Q,EL_x];
end
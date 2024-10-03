addpath ../
addpath ../Motion
clearvars

eul_sol_name = 'eul_trajectory_optimization_rigged.mat';
quat_sol_name = 'quat_trajectory_optimization_rigged.mat';
ik_sol_name = 'mot_struct_rigged.mat';

quat_struct = load(quat_sol_name);
eul_struct = load(eul_sol_name);
ik_struct = load(ik_sol_name);

quat_sol = quat2eul_motion(quat_struct.timeseries_coords.data);
eul_sol = eul_struct.timeseries_coords.data;
ik_sol = ik_struct.coords_struct.mot_euler;
time_ik = ik_struct.coords_struct.tout;
time = eul_struct.timeseries_coords.tout;

for i = 1:length(eul_sol(1,:))
    p = polyfit(time_ik,ik_sol(:,i),6);
    ik_polyfitted = polyval(p,time);

    figure
    plot(time,eul_sol(:,i),time,quat_sol(:,i),time,ik_polyfitted)
    legend('eul','quat','ik')
end


addpath Functions\

python_sol_quat = load('quat_trajectory_optimization.mat');
data2mot(python_sol_quat.timeseries_coords,'quat_trajectory_optimization.mot', 'quaternion', 'struct')

python_sol_eul = load('eul_trajectory_optimization.mat');
data2mot(python_sol_eul.timeseries_coords,'eul_trajectory_optimization.mot', 'euler', 'struct')


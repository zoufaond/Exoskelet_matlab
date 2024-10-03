addpath Functions\

% python_sol_quat = load('quat_steering_eulIC.mat');
% data2mot(python_sol_quat.timeseries_coords,'quat_steering_eulIC.mot', 'quaternion', 'struct')

python_sol_eul = load('eul_steering_eulIC.mat');
data2mot(python_sol_eul.timeseries_coords,'eul_steering_eulIC.mot', 'euler', 'struct')


clearvars

addpath Geometry\
% addpath EOMs_quat\
% addpath EOMs_eul\
addpath equations_of_motion\quaternion\
addpath equations_of_motion\euler\
addpath Muscle_full_model\euler\
addpath Muscle_full_model\quaternion\
addpath Functions\
addpath Motion\

data = load('data_model.mat');
% opensim_model = load('das3_quat.mat');
% motion = load('mot_struct.mat');
model = data.params.model;
params = data.params;

initEul = params.InitPosOptEul.initCondEul;
initQuat = params.InitPosOptEul.initCondQuat;
% initEul(8) = 0; % gimbal lock
% initQuat = motion.coords_struct.mot_quaternion_mod(1,:)';

% muscles = opensim_model.model_full_quat.muscles;
% for i=1:length(muscles)
%     fmax(i,1) = muscles{i}.fmax;%*muscles_compensations(i);
%     lceopt(i,1) = muscles{i}.lceopt;
%     lslack(i,1) = muscles{i}.lslack;
% end




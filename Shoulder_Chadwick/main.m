% clearvars

addpath Geometry\
addpath EOMs_quat\
addpath EOMs_eul\
addpath Muscle_full_model\euler\
addpath Muscle_full_model\quaternion\
addpath Functions\
data = load('data_model.mat');
opensim_model = load('das3_full_quat.mat');
model = data.params.model;
params = data.params;

initEul = data.params.initCond_Eul;
initQuat = params.initCond_Q;

muscles = opensim_model.model_full_quat.muscles;

for i=1:length(muscles)
    fmax(i,1) = muscles{i}.fmax;%*muscles_compensations(i);
    lceopt(i,1) = muscles{i}.lceopt;
    lslack(i,1) = muscles{i}.lslack;
end




% clearvars

addpath Geometry\
addpath EOMs_quat\
addpath EOMs_eul\
addpath Muscle_simplified_model\euler\
addpath Muscle_simplified_model\quaternion\
data = load('data_model_python2.mat');
model = data.params.model;
params = data.params;

initEul = data.params.InitPosOptEul.initCondEul;
initEul(8) = initEul(8) + 0.2;
initQuat = params.InitPosOptQuat.initCondQuat;
initSimScape = initEul;
% initQuat = model.initCond_Q;

% %[initCondEul;zeros(10,1)]
% simpl = load('das3_simplified_eul.mat');
% full = load('das3_full_eul.mat');
% 
% mussimp = simpl.model_simplified_eul.muscles;
% % musfull = full.model_full_eul.muscles;
% 
% for i = 1:length(mussimp)
%     simname{i,1} = mussimp{i}.name;
% end
% 
% for i = 1:length(musfull)
%     fullname{i,1} = musfull{i}.name;
% end


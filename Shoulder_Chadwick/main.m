clearvars

addpath Geometry\
addpath EOMs_quat\
addpath EOMs_eul\
addpath Muscle_simplified_model\euler\
addpath Muscle_simplified_model\quaternion\
data = load('data_model_mod.mat');
model = data.model;

initEul = model.InitPosOptEul.initCondEul;
initQuat = model.InitPosOptQuat.initCondQuat;
initSimScape = initEul;
% initQuat = model.initCond_Q;
act = zeros(36,1);
% act(6) = 0.2;
% act(20:21) = 0.2;
% act(29) = 0.2;

%[initCondEul;zeros(10,1)]
simpl = load('das3_simplified_eul.mat');
full = load('das3_full_eul.mat');

mussimp = simpl.model_simplified_eul.muscles;
% musfull = full.model_full_eul.muscles;

for i = 1:length(mussimp)
    simname{i,1} = mussimp{i}.name;
end
% 
% for i = 1:length(musfull)
%     fullname{i,1} = musfull{i}.name;
% end


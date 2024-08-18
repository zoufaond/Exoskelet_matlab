addpath euler\
addpath quaternion\
addpath ..\

% this runs functions that creates external moments from muscle forces

genEq = 1;
modelEul = load("das3_simplified_eul.mat");
modelQuat = load("das3_simplified_quat.mat");
[~,~,~] = muscle_derivation_euler(modelEul.model_simplified_eul,genEq);
% [~,~,~] = muscle_derivation_quat(modelQuat.model_simplified_quat,genEq);
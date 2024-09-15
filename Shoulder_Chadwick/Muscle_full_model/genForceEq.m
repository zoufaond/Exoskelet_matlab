addpath euler\
addpath quaternion\
addpath ..\

genEq = 1;
modelEul = load("das3_full_eul_03.mat");
modelQuat = load("das3_full_quat_03.mat");
% [~,~,~] = muscle_derivation_euler(modelEul.model_full_eul,genEq);
[~,~,~] = muscle_derivation_quat(modelQuat.model_full_quat,genEq);
addpath euler\
addpath quaternion\
addpath ..\

genEq = 1;
modelEul = load("das3_simplified.mat");
modelQuat = load("das3_simplified_quaternion.mat");
[~,~,~] = muscle_derivation_euler(modelEul.model_simpl,genEq);
% [~,~,~] = muscle_derivation_quat(modelQuat.model_simpl_quat,genEq);
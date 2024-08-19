addpath euler\
addpath quaternion\
addpath ..\
addpath ..\DAS3-simplified\model\polys_mod_full\
addpath ..\DAS3-simplified\model

genEq = 1;
modelEul = load("das3_full_eul.mat");
modelQuat = load("das3_full_quat.mat");
[~,~,~] = muscle_derivation_euler(modelEul.model_full_eul,genEq);
[~,~,~] = muscle_derivation_quat(modelQuat.model_full_quat,genEq);
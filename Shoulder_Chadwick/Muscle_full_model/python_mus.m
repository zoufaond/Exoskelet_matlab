addpath euler\
addpath quaternion\
addpath ..\
% clearvars

data = load('data_model.mat');
model = data.params.model;
params = data.params;

fmax = params.InitPosOptQuat.fmax;
lceopt = params.InitPosOptQuat.lceopt;
lslack = params.InitPosOptQuat.lslack;

acts = ones(137,1)*0.5;
q = [0,0,0,0,linspace(0.1,0.4,13),0]';
res = TE_full_quat(0,q,acts,fmax,lceopt,lslack)

% clear all
addpath Geometry\
addpath EOMs_quat\
addpath EOMs_eul\
data = load('data_model.mat');
% load('data_model_eul.mat')
initCondEul = [data.model.initCond_eul.SC;
               data.model.initCond_eul.AC;
               data.model.initCond_eul.GH;
               data.model.initCond_eul.EL];

model = data.model;



activ = zeros(35,1);



addpath('..\equations_of_motion\')
addpath('..\Muscle_simplified_model')
addpath('..\Functions\')
addpath('..\')
load('data_model.mat')
load("das3_simplified.mat")

load('info_opt_control_new.mat')

matrix = Humerus_to_Thorax_matrix(info.Xopt);
angles_res = rotxyz_sym(info.Xopt');

data_to_mot([info.Topt(1:end-1),info.Xopt(1:end-1,1:10)]','opt_contr2.sto')
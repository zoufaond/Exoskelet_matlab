% clear all
addpath('Functions\')
addpath('equations_of_motion\')
addpath('Muscle_simplified_model\')
load('data_model.mat')
load("das3_simplified.mat")
SC_yzx = [-21.784 6.303 0]*pi/180;
AC_yzx = [46.295 4.899 -1.180]*pi/180;
GH_yzy = [0 5 0]*pi/180;
EL_x = 10*pi/180;
PS_y = 5*pi/180;
initCond = [SC_yzx,AC_yzx,GH_yzy,EL_x,zeros(1,10)]';
muscles = model_simpl.muscles;
for i=1:length(muscles)
    name_vec{i,1} = muscles{i}.name;
    fmax_vec(i,1) = muscles{i}.fmax;
    lceopt_vec(i,1) = muscles{i}.lceopt;
    lslack_vec(i,1) = muscles{i}.lslack;
end
% % % % 
% fmax_optim = model.q_fmax_lceopt_InOut.fmax_optim;
% lceopt_optim = model.q_fmax_lceopt_InOut.lceopt_optim;
% lslack_optim = model.q_fmax_lceopt3.lslack_optim;
% initCond_optim = model.q_fmax_lceopt3.initCond_optim;
activ = zeros(35,1);
% activ(10:12) = 0.1;

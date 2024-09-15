clearvars
addpath EOMs_eul\
addpath Muscle_full_model\euler\
addpath Functions\
addpath Motion\

opensim_model = load('das3_full_eul_03.mat');
data = load('data_model.mat');
motion = load('mot_struct.mat');
initCond = motion.coords_struct.mot_euler(1,:)';
% 
muscles = opensim_model.model_full_eul.muscles;
% 
for i=1:length(muscles)
    fmax_vec(i,1) = muscles{i}.fmax;
    lceopt_vec(i,1) = muscles{i}.lceopt;
    lslack_vec(i,1) = muscles{i}.lslack;
end
%%
t = 0;
qn = 10;
nmus = length(fmax_vec);

fun = @(x) sum(mus_forces_eul(t,x(1:qn,1), ...
    zeros(nmus,1),x(qn+1:qn+nmus,1),x(qn+nmus+1:(end-1),1),lslack_vec).^2);
  
x0 = [initCond(1:qn);fmax_vec;lceopt_vec;1]; %;
A = [];
b = [];
Aeq = [];
beq = [];
AC_bndrs = [ones(3,1)*0.005;ones(3,1)*0.005;ones(3,1)*0.005;0.005];
lb = [initCond(1:10)-AC_bndrs;fmax_vec-fmax_vec*0.1;lceopt_vec-lceopt_vec*0.1;0.5];
ub = [initCond(1:10)+AC_bndrs;fmax_vec+fmax_vec*0.1;lceopt_vec+lceopt_vec*0.1;1.5];
nonlcon = @(x) moment_equilibrium(t,x(1:qn,1),zeros(nmus,1),x(qn+1:qn+nmus,1),x(qn+nmus+1:(end-1),1),lslack_vec,data.params.model,x(end));

options = optimoptions(@fmincon,'Display','iter','MaxFunEval',1e7,'algorithm','interior-point','MaxIter',10000);%,'MaxFunEval',1e7,'TolFun',1e-9,'MaxIter',1e6,'algorithm','interior-point','TolCon',1e-8,'TolX',1e-12);
x_new = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
%%
fmax_optim = x_new(11:11+nmus-1);
lceopt_optim = x_new(11+nmus:end-1);
lslack_optim = lslack_vec;
initCond_optim = [x_new(1:10);zeros(10,1)];
data.params.model.first_elips_scale = x_new(end);
data.params.model.second_elips_scale = x_new(end)+0.2;
data.params.InitPosOptEul.fmax = fmax_optim;
data.params.InitPosOptEul.lceopt = lceopt_optim;
data.params.InitPosOptEul.lslack = lslack_optim;
data.params.InitPosOptEul.initCondEul = initCond_optim;
eul2quatConv = eul2quat_motion(initCond_optim');
data.params.InitPosOptEul.initCondQuat = [eul2quatConv';zeros(10,1)];
% save('data_model_mod.mat','model')
%%
get_conoid_length(initCond_optim(4:6),data.params.model)
%%
params = data.params;
save('data_model.mat','params')
%%
eul2quatConv = eul2quat_motion(data.params.InitPosOptEul.initCondEul');
data.params.InitPosOptEul.initCondQuat = [eul2quatConv';zeros(10,1)];
params = data.params;
save('data_model.mat','params')

%%
function [c,ceq] = moment_equilibrium(t,q,act,fmax_vec,lceopt_vec,lslack_vec,model,first_elips_scale)
    x = q(1:10);
    FO = fo_EUL(t,x,zeros(10,1),act,model,first_elips_scale);
    FE_muscles = TE_eul(t,x,act,fmax_vec,lceopt_vec,lslack_vec);
    ceq = FO+[zeros(10,1);FE_muscles];
    c = [];
end
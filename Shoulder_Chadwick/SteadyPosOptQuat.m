addpath EOMS_quat
addpath Muscle_full_model\quaternion\
addpath Functions\
addpath Motion\

%%
clearvars

opensim_model = load('das3_full_quat.mat');
data = load('data_model.mat');
motion = load('mot_struct.mat');
initCond = motion.coords_struct.mot_quaternion(1,:)';
%%
muscles = opensim_model.model_full_quat.muscles;

for i=1:length(muscles)
    fmax_vec(i,1) = muscles{i}.fmax;%*muscles_compensations(i);
    lceopt_vec(i,1) = muscles{i}.lceopt;
    lslack_vec(i,1) = muscles{i}.lslack;
end
%%

t = 0;
qn = 13;
nmus = length(fmax_vec);

fun = @(x) 1 * sum(mus_forces_quat(t,x(1:qn,1), ...
    zeros(nmus,1),x(qn+1:qn+nmus,1),x(qn+nmus+1:(end-1),1),lslack_vec).^2);
  
x0 = [initCond(1:qn);fmax_vec;lceopt_vec;1]; %;
A = [];
b = [];
Aeq = [];
beq = [];
AC_bndrs = [ones(4,1)*0.01;ones(4,1)*0.01;ones(4,1)*0.01;0.005];
lb = [initCond(1:qn)-AC_bndrs;fmax_vec-fmax_vec*0.2;lceopt_vec-lceopt_vec*0.15;0.5];
ub = [initCond(1:qn)+AC_bndrs;fmax_vec+fmax_vec*0.2;lceopt_vec+lceopt_vec*0.15;1.5];
% nonlcon =@(x) [0,0];
nonlcon = @(x) moment_equilibrium(t,x(1:qn,1),zeros(nmus,1),x(qn+1:qn+nmus,1),x(qn+nmus+1:(end-1),1),lslack_vec,data.params.model,x(end));

options = optimoptions(@fmincon,'Display','iter','MaxFunEval',1e7,'algorithm','interior-point','MaxIter',10000);%,'MaxFunEval',1e7,'TolFun',1e-9,'MaxIter',1e6,'algorithm','interior-point','TolCon',1e-8,'TolX',1e-12);
x_new = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
 %%
% % % options.Algorithm = 'interior-point';
% options.Algorithm = 'interior-point';
% options.MaxIterations = 500;
% x_new = fmincon(fun,x_new,A,b,Aeq,beq,lb,ub,nonlcon,options);
%%
fmax_optim = x_new(qn+1:qn+1+nmus-1);
lceopt_optim = x_new(qn+1+nmus:end);
lslack_optim = lslack_vec;
initCond_optim = [x_new(1:qn);zeros(10,1)];
data.params.InitPosOptQuat.fmax = fmax_optim;
data.params.InitPosOptQuat.lceopt = lceopt_optim;
data.params.InitPosOptQuat.lslack = lslack_optim;
data.params.InitPosOptQuat.initCondQuat = initCond_optim;
quat2eulConv = quat2eul_motion(initCond_optim');
data.params.InitPosOptQuat.initCondEul = [quat2eulConv';zeros(10,1)];
% 
params = data.params;
save('data_model.mat','params')
%%
function [c,ceq] = moment_equilibrium(t,q,act,fmax_vec,lceopt_vec,lslack_vec,model,opt_var)

    FO = fo_quat(t,q(1:13),zeros(10,1),model,[opt_var]);
    FE_muscles = TE_quat(t,q,act,fmax_vec,lceopt_vec,lslack_vec);
    ceq = [FO+[zeros(13,1);FE_muscles];sum(q(1:4).^2)-1;sum(q(5:8).^2)-1;sum(q(9:12).^2)-1];
    c = [];
end
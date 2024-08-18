clearvars
addpath EOMS_quat
addpath Muscle_simplified_model\quaternion\
addpath Functions\
% addpath('Muscle_full_model\')

model_quat = load('das3_simplified_quat.mat');
data = load('data_model_python.mat');
initCond = data.params.initCond_Q;
%%
muscles = model_quat.model_simplified_quat.muscles;

muscles_compensations = [11/5,11/5,2,2,4,5,4,4,4,11/3,11/3,11/3,4,3,6,3,4,4,11,1,1,4,3,3,3,3,2,5,7,3,1,1,5,3,5,5];

for i=1:length(muscles)
    fmax_vec(i,1) = muscles{i}.fmax*muscles_compensations(i);
    lceopt_vec(i,1) = muscles{i}.lceopt;
    lslack_vec(i,1) = muscles{i}.lslack;
end
%%
% mus_forces_simp_ulna_quat(t,[zeros(4,1);initCond(1:qn);0],zeros(nmus,1),fmax_vec,lceopt_vec,lslack_vec)
%%
t = 0;
qn = 13;
nmus = length(fmax_vec);


fun = @(x) sum(mus_forces_simplified_quat(t,[zeros(4,1);x(1:qn,1);0], ...
    zeros(nmus,1),x(qn+1:qn+nmus,1),x(qn+nmus+1:end,1),lslack_vec).^2);
  
x0 = [initCond(1:qn);fmax_vec;lceopt_vec]; %;
A = [];
b = [];
Aeq = [];
beq = [];
AC_bndrs = [ones(4,1)*0.05;ones(4,1)*0.05;ones(4,1)*0.005;0.05];
lb = [initCond(1:qn)-AC_bndrs;fmax_vec-fmax_vec*0.2;lceopt_vec-lceopt_vec*0.25];
ub = [initCond(1:qn)+AC_bndrs;fmax_vec+fmax_vec*0.2;lceopt_vec+lceopt_vec*0.1];
% nonlcon =@(x) [0,0];
nonlcon = @(x) moment_equilibrium(t,x(1:qn,1),zeros(nmus,1),x(qn+1:qn+nmus,1),x(qn+nmus+1:end,1),lslack_vec,data.params.model);

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
data.params.InitPosOptQuat.initCondEul = [quat2eulConv(1:6),data.params.initCond_Eul(7:10)',zeros(1,10)]';

params = data.params;
save('data_model_python.mat','params')
%%
function [c,ceq] = moment_equilibrium(t,q,act,fmax_vec,lceopt_vec,lslack_vec,model)

    FO = fo_Q(t,q(1:13),zeros(10,1),model);
    % FE_muscles = f_muscles_full_ulna(t,x(1:10),act,fmax_vec,lceopt_vec,lslack_vec);
    FE_muscles = TE_simplified_quat(t,[zeros(4,1);q;0],act,fmax_vec,lceopt_vec,lslack_vec);
    ceq = [FO+FE_muscles;sum(q(1:4).^2)-1;sum(q(5:8).^2)-1;sum(q(9:12).^2)-1];
    c = [];
end
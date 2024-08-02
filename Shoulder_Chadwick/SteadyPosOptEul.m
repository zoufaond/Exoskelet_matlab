clearvars
addpath EOMs_eul\
addpath Muscle_simplified_model\euler\
data = load('data_model_mod.mat');

initCond = data.model.InitPosOptQuat.initCondEul;
opensim_model = load('das3_simplified_eul.mat');
% %%
% muscles = opensim_model.model_simplified_eul.muscles;
% 
% muscles_compensations = [11/5,11/5,2,2,4,5,4,4,4,11/3,11/3,11/3,4,3,6,3,4,4,11,1,1,4,3,3,3,3,2,5,7,3,1,1,5,3,5,5];
% 
% for i=1:length(muscles)
%     fmax_vec(i,1) = muscles{i}.fmax*muscles_compensations(i);
%     lceopt_vec(i,1) = muscles{i}.lceopt;
%     lslack_vec(i,1) = muscles{i}.lslack;
% end
fmax_vec = data.model.InitPosOptQuat.fmax;
lceopt_vec = data.model.InitPosOptQuat.lceopt;
lslack_vec = data.model.InitPosOptQuat.lslack;
%%
t = 0;
qn = 10;
nmus = length(fmax_vec);


fun = @(x) sum(mus_forces_simplified_eul(t,[zeros(3,1);x(1:qn,1);0], ...
    zeros(nmus,1),x(qn+1:qn+nmus,1),x(qn+nmus+1:end,1),lslack_vec).^2);

% fun = @(x) sum(moment_equilibrium(t,x(1:qn,1),zeros(nmus,1),x(qn+1:qn+nmus,1),x(qn+nmus+1:end,1),lslack_vec,model).^2);
    
x0 = [initCond(1:10);fmax_vec;lceopt_vec]; %;
A = [];
b = [];
Aeq = [];
beq = [];
AC_bndrs = [ones(3,1)*0.2;ones(3,1)*0.2;ones(3,1)*0.2;0.01];
lb = [initCond(1:10)-AC_bndrs;fmax_vec-fmax_vec*0.3;lceopt_vec-lceopt_vec*0.4];
ub = [initCond(1:10)+AC_bndrs;fmax_vec+fmax_vec*0.3;lceopt_vec+lceopt_vec*0.4];
nonlcon = @(x) moment_equilibrium(t,x(1:qn,1),zeros(nmus,1),x(qn+1:qn+nmus,1),x(qn+nmus+1:end,1),lslack_vec,data.model);

options = optimoptions(@fmincon,'Display','iter','MaxFunEval',1e7,'algorithm','interior-point','MaxIter',10000);%,'MaxFunEval',1e7,'TolFun',1e-9,'MaxIter',1e6,'algorithm','interior-point','TolCon',1e-8,'TolX',1e-12);
x_new = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
%  %%
% % % options.Algorithm = 'interior-point';
% options.Algorithm = 'interior-point';
% options.MaxIterations = 1000;
% x_new = fmincon(fun,x_new,A,b,Aeq,beq,lb,ub,nonlcon,options);
%%
fmax_optim = x_new(11:11+nmus-1);
lceopt_optim = x_new(11+nmus:end);
lslack_optim = lslack_vec;
initCond_optim = [x_new(1:10);zeros(10,1)];
data.model.InitPosOptEul.fmax = fmax_optim;
data.model.InitPosOptEul.lceopt = lceopt_optim;
data.model.InitPosOptEul.lslack = lslack_optim;
data.model.InitPosOptEul.initCondEul = initCond_optim;

model = data.model;
save('data_model_mod.mat','model')
%%
% model = data.model;
% save('data_model.mat','model')
%%
function [c,ceq] = moment_equilibrium(t,q,act,fmax_vec,lceopt_vec,lslack_vec,model)
    x = q(1:10);
    FO = fo_EUL(t,x,zeros(10,1),model);
    FE_muscles = TE_simplified_eul(t,[zeros(3,1);x(1:10);0],act,fmax_vec,lceopt_vec,lslack_vec);
    ceq = FO+FE_muscles;
    c = [];
end
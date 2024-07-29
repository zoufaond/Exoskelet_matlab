addpath('..\')
clearvars
clc
model_eul = load('data_model_eul.mat');
model_quat = load('data_model_quat.mat');
x0 = model_eul.model.q_fmax_lceopt_InOut2.initCond_optim(1:10);
% ID_fun([x0;zeros(10,1);zeros(10,1)],zeros(10,1),zeros(10,1),zeros(35,1),model)

scale = 0.2;
n = 2;
t_end = 10;
[EulTraj,QuatTraj] = create_abduction_traj(x0,n,t_end,scale);
%%
x0 = ones(35,1)*0.0;
for i=1:n
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    lb = zeros(35,1);
    ub = ones(35,1)*100;
    fun = @(x) sum(x.^2);
    options = optimoptions(@fmincon,'Display','iter','MaxFunEval',1e7,'algorithm','SQP','MaxIter',10000);
    nonlcon = @(x) ID_Quat_cceq(QuatTraj.quat(:,i),QuatTraj.dquat(:,i),QuatTraj.omega(:,i),QuatTraj.domega(:,i),x(:,1),model_quat.model);
    x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
    actQuat(:,i) = x;
    x0 = x;
end

%%
x0 = ones(35,1)*0.0;
for i=1:n
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    lb = zeros(35,1);
    ub = ones(35,1)*1;
    fun = @(x) sum(x.^2);
    options = optimoptions(@fmincon,'Display','iter','MaxFunEval',1e7,'algorithm','SQP','MaxIter',10000);
    nonlcon = @(x) ID_Eul_cceq(EulTraj.q(:,i),EulTraj.dq(:,i),EulTraj.ddq(:,i),x(:,1),model_eul.model);
    x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
    actEul(:,i) = x;
    x0 = x;
end
% %%
% InitGuess.T = linspace(0,t_end,n);
% InitGuess.MV0 = x;
% InitGuess.X0 = [q;dq];
% save('InitGuess.mat','InitGuess')
% clearvars
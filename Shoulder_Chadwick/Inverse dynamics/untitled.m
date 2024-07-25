addpath('..\equations_of_motion\')
addpath('..\Muscle_simplified_model')
addpath('..\Functions\')
addpath('..\')
load('data_model.mat')
x0 = model.q_fmax_lceopt_InOut2.initCond_optim(1:10);
% ID_fun([x0;zeros(10,1);zeros(10,1)],zeros(10,1),zeros(10,1),zeros(35,1),model)

scale = 0.5;
n = 8;
t_end = 1;
[traj,q,dq,ddq] = create_abduction_traj(x0,n,t_end,scale);
%%
for i=1:n
x0 = ones(35,1)*0; %;
A = [];
b = [];
Aeq = [];
beq = [];
lb = zeros(35,1);
ub = ones(35,1)*100;
fun = @(x) sum(x.^2);
options = optimoptions(@fmincon,'MaxFunEval',1e7,'algorithm','interior-point','MaxIter',10000);
nonlcon = @(x) ID_fun(q(:,i),dq(:,i),ddq(:,i),x(:,1),model);
x(:,i) = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
end
%%
InitGuess.T = linspace(0,t_end,n);
InitGuess.MV0 = x;
InitGuess.X0 = [q;dq];
save('InitGuess.mat','InitGuess')
clearvars
% addpath('Functions\')
addpath equations_of_motion\
% addpath('Muscle_full_model\')
addpath('Muscle_simplified_model\')
load('data_model.mat')

SC_yzx = [-21.784 6.303 0]*pi/180;
AC_yzx = [46.295 4.899 -1.180]*pi/180;
GH_yzy = [0 5 0]*pi/180;
EL_x = 10*pi/180;
PS_y = 5*pi/180;
initCond = [SC_yzx,AC_yzx,GH_yzy,EL_x,zeros(1,10)]';
load('das3_simplified.mat')
%%
muscles = model_simpl.muscles;
for i=1:length(muscles)
    fmax_vec(i,1) = muscles{i}.fmax*4;
    lceopt_vec(i,1) = muscles{i}.lceopt;
    lslack_vec(i,1) = muscles{i}.lslack;
end
%%
t = 0;
qn = 10;
nmus = length(fmax_vec);


fun = @(x) sum(mus_forces_ulna(t,x(1:qn,1), ...
    zeros(nmus,1),x(qn+1:qn+nmus,1),x(qn+nmus+1:end,1),lslack_vec).^2);

% fun = @(x) sum(moment_equilibrium(t,x(1:qn,1),zeros(nmus,1),x(qn+1:qn+nmus,1),x(qn+nmus+1:end,1),lslack_vec,model).^2);
    
x0 = [initCond(1:10);fmax_vec;lceopt_vec]; %;
A = [];
b = [];
Aeq = [];
beq = [];
AC_bndrs = [ones(3,1)*0.01;ones(3,1)*0.15;ones(4,1)*0.005];
lb = [initCond(1:10)-AC_bndrs;fmax_vec-fmax_vec*0.5;lceopt_vec-lceopt_vec*0.7];
ub = [initCond(1:10)+AC_bndrs;fmax_vec+fmax_vec*10000;lceopt_vec+lceopt_vec*0.05];
% nonlcon =@(x) [0,0];
nonlcon = @(x) moment_equilibrium(t,x(1:qn,1),zeros(nmus,1),x(qn+1:qn+nmus,1),x(qn+nmus+1:end,1),lslack_vec,model);

options = optimoptions(@fmincon,'Display','iter','MaxFunEval',1e7,'algorithm','interior-point','MaxIter',10000);%,'MaxFunEval',1e7,'TolFun',1e-9,'MaxIter',1e6,'algorithm','interior-point','TolCon',1e-8,'TolX',1e-12);
x_new = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
 %%
% % options.Algorithm = 'interior-point';
% options.Algorithm = 'sqp';
% options.MaxIterations = 600;
% x_new = fmincon(fun,x_new,A,b,Aeq,beq,lb,ub,nonlcon,options);
%%
fmax_optim = x_new(11:11+nmus-1);
lceopt_optim = x_new(11+nmus:end);
lslack_optim = lslack_vec;
initCond_optim = [x_new(1:10);zeros(10,1)];
model.q_fmax_lceopt_InOut2.fmax_optim = fmax_optim;
model.q_fmax_lceopt_InOut2.lceopt_optim = lceopt_optim;
model.q_fmax_lceopt_InOut2.lslack_optim = lslack_optim;
model.q_fmax_lceopt_InOut2.initCond_optim = initCond_optim;
save('data_model.mat','model')
%%
function [c,ceq] = moment_equilibrium(t,q,act,fmax_vec,lceopt_vec,lslack_vec,model)
    x = [q(1:10);zeros(10,1)];
    FO = fo_py_ulna_InOutCont(t,x,model);
    % FE_muscles = f_muscles_full_ulna(t,x(1:10),act,fmax_vec,lceopt_vec,lslack_vec);
    FE_muscles = f_muscles_ulna(t,x(1:10),act,fmax_vec,lceopt_vec,lslack_vec);
    ceq = FO+FE_muscles;
    c = [];
end
% clc
% clear all
addpath('..\equations_of_motion\')
addpath('..\Muscle_simplified_model')
addpath('..\Functions\')
addpath('..\Inverse Dynamics\')
addpath('..\')
load('data_model.mat')
load("das3_simplified.mat")
load('InitGuess')
%%
nx = 20;
ny = 20;
nu = 35;
nlobj = nlmpc(nx,ny,nu);
Ts = 1/3;
p_hor = 3;
c_hor = 3;
sim_time = Ts*p_hor;
nlobj.Ts = Ts;
x0 = model.q_fmax_lceopt_InOut2.initCond_optim;
for i=1:nu
nlobj.MV(i).Min = 0;
nlobj.MV(i).Max = 1;
end

for i=1:10
    nlobj.States(i).Min = -1;
    nlobj.States(i).Max = 1;
end

for i=11:20
    nlobj.States(i).Min = -5;
    nlobj.States(i).Max = 5;
end


[~,q_traj,~,~] = create_abduction_traj(x0,p_hor,p_hor*Ts,0.2);
%%
nlobj.PredictionHorizon = p_hor;
nlobj.ControlHorizon = c_hor;
nlobj.Model.StateFcn = "nlmpc_fun";
% +(sum((X(10:p_hor,1)-0).^2)+sum((X(10:p_hor,2)-0.3).^2))*100

nlobj.Model.NumberOfParameters = 1;
nlobj.Optimization.CustomCostFcn = @(X,U,e,data,model) sum(sum((X(2:end,1:10)-q_traj').^2));% + 1e-4*sum(sum(U(1:end,:).^2));
% nlobj.Optimization.CustomCostFcn = @(X,U,e,data,model) sum(sum((rotxyz_sym(X(end,:)')-traj).^2)) + sum(sum(U(1:p_hor,:).^2));
% nlobj.Optimization.CustomEqConFcn = @(X,U,data,model) ;
nlobj.Optimization.ReplaceStandardCost = true;
nlobj.Optimization.SolverOptions.Display = "iter-detailed";
nlobj.Optimization.SolverOptions.MaxIterations = 1e4;
nlobj.Optimization.SolverOptions.StepTolerance = 1e-20;
% nlobj.Optimization.SolverOptions.OptimalityTolerance = 1e-5;
% nlobj.Optimization.SolverOptions.ConstraintTolerance = 1e-5;
% nlobj.Optimization.SolverOptions.FunctionTolerance = 1e-6;
nlobj.Optimization.SolverOptions.MaxFunctionEvaluations = 1e7;
nlobj.Optimization.SolverOptions.Algorithm = "interior-point";
nlobj.Optimization.SolverOptions.UseParallel = true;

nloptions = nlmpcmoveopt;
nloptions.Parameters = {model};

%%
nloptions.X0 = [InitGuess.X0(:,1),InitGuess.X0]';
nloptions.MV0 =[InitGuess.MV0(:,1),InitGuess.MV0]';
u0 = InitGuess.MV0(:,1)*0;
validateFcns(nlobj,x0,u0,[],{model});
%%

[mv,simdata,info] = nlmpcmove(nlobj,x0,u0,[],[],nloptions);

%%
% data_to_mot([info.Topt,info.Xopt(:,1:10)]','opt_contr.sto')
% save('info_opt_control','info')
% inputData6 = timeseries(info.MVopt(1:end-pre_end,6),info.Topt(1:end-pre_end));
% save("inputData1.mat","inputData1","-v7.3");
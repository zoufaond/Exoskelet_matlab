% clc
% clear all
addpath('..\equations_of_motion\')
addpath('..\Muscle_simplified_model')
addpath('..\Functions\')
addpath('..\')
load('data_model.mat')
load("das3_simplified.mat")
%%
nx = 20;
ny = 20;
nu = 35;
nlobj = nlmpc(nx,ny,nu);
Ts = 0.1;
p_hor = 5;
c_hor = 5;
sim_time = Ts*p_hor;
nlobj.Ts = Ts;
x0 = model.q_fmax_lceopt_InOut.initCond_optim;
for i=1:nu
nlobj.MV(i).Min = 0;
nlobj.MV(i).Max = 1e4;
end

traj = create_abduction_traj(model.q_fmax_lceopt_InOut2.initCond_optim,p_hor);
%%
nlobj.PredictionHorizon = p_hor;
nlobj.ControlHorizon = c_hor;
nlobj.Model.StateFcn = "nlmpc_fun";
% +(sum((X(10:p_hor,1)-0).^2)+sum((X(10:p_hor,2)-0.3).^2))*100

nlobj.Model.NumberOfParameters = 1;
% nlobj.Optimization.CustomCostFcn = @(X,U,e,data,model) sum(sum((Humerus_to_Thorax_matrix(rotxyz(X(end,:)))-traj(:,end)').^2)) + sum(sum(U(1:p_hor,:).^2));
nlobj.Optimization.CustomCostFcn = @(X,U,e,data,model) sum( (rotxyz(Humerus_to_Thorax_matrix(X(end,:)))-traj(:,end)').^2 );
% nlobj.Optimization.CustomEqConFcn = @(X,U,data,model) ;
nlobj.Optimization.ReplaceStandardCost = true;
nlobj.Optimization.SolverOptions.Display = "iter-detailed";
nlobj.Optimization.SolverOptions.MaxIterations = 1e4;
% nlobj.Optimization.SolverOptions.StepTolerance = 1e-7;
% nlobj.Optimization.SolverOptions.OptimalityTolerance = 1e-5;
% nlobj.Optimization.SolverOptions.ConstraintTolerance = 1e-5;
% nlobj.Optimization.SolverOptions.FunctionTolerance = 1e-6;
nlobj.Optimization.SolverOptions.MaxFunctionEvaluations = 1e7;
nlobj.Optimization.SolverOptions.Algorithm = "interior-point";
nlobj.Optimization.SolverOptions.UseParallel = true;

nloptions = nlmpcmoveopt;
nloptions.Parameters = {model};
%%
load("info_opt_control.mat")
Xopt0 = opt_control_1.Xopt;
MVopt0 = opt_control_1.MVopt*2;
% abduct = [zeros(8,1);0.1;0;zeros(10,1)];
% nloptions.X0 = [x0,x0+abduct,x0+abduct*2]';
% musclex0 = zeros(35,1);
% musclex0(10:12) = 1;
%%
nloptions.X0 = Xopt0;
nloptions.MV0 = MVopt0;
u0 = ones(35,1)*0;
validateFcns(nlobj,x0,u0,[],{model});
%%

[mv,simdata,info] = nlmpcmove(nlobj,x0,u0,[],[],nloptions);

%%
% data_to_mot([info.Topt,info.Xopt(:,1:10)]','opt_contr.sto')
% save('info_opt_control','info')
% inputData6 = timeseries(info.MVopt(1:end-pre_end,6),info.Topt(1:end-pre_end));
% save("inputData1.mat","inputData1","-v7.3");
addpath('..\equations_of_motion\')
addpath('..\Muscle_simplified_model')
addpath('..\Functions\')
addpath('..\Inverse dynamics\')
addpath('..\..\..\OptimTraj\')
load('data_model.mat')
load("das3_simplified.mat")
load('InitGuess.mat')

t0 = 0;
tF = 1;
N_grid = 8;
x0 = model.q_fmax_lceopt_InOut2.initCond_optim;
scale = 0.5;
% [~,q_traj,~,~] = create_abduction_traj(x0,N_grid,0.1,0.1);
%%

problem.func.dynamics = @(t,x,u)( OptimTraj_fun(x,u,model) );
problem.func.pathObj = @(t,x,u) sum((reference_trajectory(t,tF,x0,scale)-x(1:10,:)).^2);% + 1e-10*(sum(u.^2));

problem.bounds.state.low = -ones(20,1)*10;
problem.bounds.state.upp =  ones(20,1)*10;

problem.bounds.initialState.low = x0;
problem.bounds.initialState.upp = x0;
problem.bounds.finalState.low = -inf(20,1);
problem.bounds.finalState.upp = inf(20,1);

problem.bounds.initialTime.low = t0;
problem.bounds.initialTime.upp = t0;
problem.bounds.finalTime.low = tF;
problem.bounds.finalTime.upp = tF;

problem.bounds.control.low = ones(35,1)*0;
problem.bounds.control.upp = ones(35,1)*2;

% problem.guess.time = [0,tF];
% problem.guess.state = [x0, x0];
% problem.guess.control = zeros(35,2);
problem.guess.time = InitGuess.T;
problem.guess.state = InitGuess.X0+rand(20,N_grid)*0.001;
problem.guess.control = InitGuess.MV0+rand(35,N_grid)*0.001;

problem.options.nlpOpt = optimset(...
    'Display','iter',...
    'MaxFunEvals',1e6,...
    'MaxIter', 2.5e4,...
    'algorithm','interior-point', ...
    'UseParallel',true);

problem.options.method = 'trapezoid';
problem.options.trapezoid.nGrid = N_grid;

soln = optimTraj(problem);
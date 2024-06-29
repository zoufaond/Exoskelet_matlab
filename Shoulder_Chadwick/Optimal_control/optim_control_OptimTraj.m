addpath ../../OptimTraj/
addpath EoMs_derivation\
addpath Muscle_full_model\

load('data_model.mat')
t0 = 0;
tF = 0.1;
initCond_optim = model.q_fmax_lceopt_full.initCond_optim;

x0 = initCond_optim;   %[q1;q2];  %initial angles   %Stable equilibrium
xF = [initCond_optim+[0;2*pi/180;zeros(18,1)]];
dx0 = zeros(10,1);   %[dq1;dq2];  %initial angle rates
dxF = zeros(10,1);

problem.func.dynamics = @(t,x,u)( OptimTraj_fun(x,u,model) );
problem.func.pathObj = @(t,x,u)( (sum(u.^2)) );

problem.bounds.state.low = -ones(20,1)*3;
problem.bounds.state.upp = ones(20,1)*3;

problem.bounds.initialState.low = -ones(20,1);
problem.bounds.initialState.upp = ones(20,1);
problem.bounds.finalState.low = -ones(20,1);
problem.bounds.finalState.upp = ones(20,1);

problem.bounds.initialTime.low = t0;
problem.bounds.initialTime.upp = t0;
problem.bounds.finalTime.low = tF;
problem.bounds.finalTime.upp = tF;

problem.bounds.control.low = zeros(138,1);
problem.bounds.control.upp = ones(138,1)*100;

problem.guess.time = [0,tF];
problem.guess.state = [x0, xF];
problem.guess.control = zeros(138,2)+0.0001;

problem.options.nlpOpt = optimset(...
    'Display','iter',...
    'MaxFunEvals',1e6,...
    'MaxIter', 2.5e4,...
    'algorithm','sqp');

problem.options.method = 'trapezoid';
problem.options.trapezoid.nGrid = 3;

soln = optimTraj(problem);
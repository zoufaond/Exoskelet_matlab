function dx = OptimTraj_fun(x,u,model)
act = u;
t = 0;
fmax_optim = model.q_fmax_lceopt_InOut2.fmax_optim;
lceopt_optim = model.q_fmax_lceopt_InOut2.lceopt_optim;
lslack_optim = model.q_fmax_lceopt_InOut2.lslack_optim;
nTime = length(x(1,:));

for i=1:nTime
    MM = mm_py_ulna_InOutCont(t,x(:,i),model);
    FO = fo_py_ulna_InOutCont(t,x(:,i),model);
    FE_muscles = f_muscles_ulna(t,x(:,i),act(:,i),fmax_optim,lceopt_optim,lslack_optim);
    dx(:,i) = MM\(FO+FE_muscles);
end
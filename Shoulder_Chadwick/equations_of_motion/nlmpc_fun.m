function dx = nlmpc_fun(x,u,model)
t = 0;
fmax_optim = model.q_fmax_lceopt_full.fmax_optim;
lceopt_optim = model.q_fmax_lceopt_full.lceopt_optim;
lslack_optim = model.q_fmax_lceopt_full.lslack_optim;
MM = mm_py_ulna(t,x,model);
FO = fo_py_ulna(t,x,model);
FE_muscles = f_muscles_full_ulna(t,x(1:10),u,fmax_optim,lceopt_optim,lslack_optim);
dx = MM\(FO+FE_muscles);
end
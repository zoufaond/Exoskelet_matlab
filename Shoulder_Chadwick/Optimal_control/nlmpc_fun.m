function dx = nlmpc_fun(x,u,model)
t = 0;
fmax_optim = model.q_fmax_lceopt_InOut.fmax_optim;
lceopt_optim = model.q_fmax_lceopt_InOut.lceopt_optim;
lslack_optim = model.q_fmax_lceopt_InOut.lslack_optim;
MM = mm_py_ulna_InOutCont(t,x,model);
FO = fo_py_ulna_InOutCont(t,x,model);
FE_muscles = f_muscles_ulna(t,x(1:10),u,fmax_optim,lceopt_optim,lslack_optim);
dx = MM\(FO+FE_muscles);
end
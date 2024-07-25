function [c,ceq] = ID_fun(x,dx,ddx,u,model)
t = 0;

fmax_optim = model.q_fmax_lceopt_InOut2.fmax_optim;
lceopt_optim = model.q_fmax_lceopt_InOut2.lceopt_optim;
lslack_optim = model.q_fmax_lceopt_InOut2.lslack_optim;
MM = mm_py_ulna_InOutCont(t,[x;dx],model);
FO = fo_py_ulna_InOutCont(t,[x;dx],model);
FE_mus = f_muscles_ulna(t,x,u,fmax_optim,lceopt_optim,lslack_optim);
ddq = [dx;ddx];
ceq = MM*ddq-FO-FE_mus;
c = [];
end
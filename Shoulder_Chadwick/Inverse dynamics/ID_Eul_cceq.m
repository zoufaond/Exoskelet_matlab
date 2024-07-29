function [c,ceq] = ID_Eul_cceq(x,dx,ddx,u,model)
addpath ..\Muscle_simplified_model\euler\
addpath ..\EOMs_eul\
t = 0;

% load optimized fmax, lceopt, lslack
fmax = model.q_fmax_lceopt_InOut2.fmax_optim;
lceopt = model.q_fmax_lceopt_InOut2.lceopt_optim;
lslack = model.q_fmax_lceopt_InOut2.lslack_optim;


MM = mm_py_ulna_InOutCont(t,[x;dx],model);
FO = fo_py_ulna_InOutCont(t,[x;dx],model);
FE_mus = FG_muscles_simp_ulna(t,[zeros(3,1);x;0],u,fmax,lceopt,lslack);

ddq = [dx;ddx];
ceq = MM*ddq-FO-FE_mus;
c = [];
end
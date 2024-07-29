function [c,ceq] = ID_Quat_cceq(q,dq,omega,domega,u,model)
addpath ..\Muscle_simplified_model\quaternion\
addpath ..\EOMs_quat\
t = 0;

% load optimized fmax, lceopt, lslack
fmax = model.muscle_params_opt.fmax;
lceopt = model.muscle_params_opt.lceopt;
lslack = model.muscle_params_opt.lslack;


MM = mm_Q(t,q,omega,model);
FO = fo_Q(t,q,omega,model);
FE_mus = TE_muscles_simp_ulna(t,[zeros(4,1);q;0],u,fmax,lceopt,lslack);
% 
ddq = [dq;domega];
ceq = MM*ddq-FO-FE_mus;
c = [];
end
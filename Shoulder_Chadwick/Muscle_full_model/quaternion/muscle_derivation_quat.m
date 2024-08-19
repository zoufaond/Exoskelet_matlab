function [jacf,mus_lenf,jacSize] = muscle_derivation_quat(model,genEq)
muscles = model.muscles;

syms t real
q = sym('q', [18,1]);
qpol = [q(2:4);q(6:8);q(10:12);q(14:16);q(17);q(18)];
qjac = [q(2:4);q(6:8);q(10:12);q(14:16);q(17);q(18)];
act = sym('act', [length(muscles) 1]);
fmax = sym('fmax', [length(muscles) 1]);
lceopt = sym('lceopt', [length(muscles) 1]);
lslack = sym('lslack', [length(muscles) 1]);
% 
for mus=1:length(muscles)
    name = muscles{mus}.name;
    npolterms = muscles{mus}.lparam_count;
    polcoeff = muscles{mus}.lcoefs;
    expon = muscles{mus}.lparams;
    musdof = muscles{mus}.dof_indeces;
    nmusdof = muscles{mus}.dof_count;
    L = 0; % Initialize the muscle length

    for i = 1:npolterms
        % Add this term's contribution to the muscle length
        term = polcoeff(i);
        for j = 1:nmusdof
            mdof = musdof(j);
            for k = 1:expon(i, j)
                term = term * qpol(mdof);
            end
        end
        L = L + term;
    end
    mus_len(mus) = L;

    if genEq == 1
        mus_force(mus) = muscle_force(act(mus),L,fmax(mus),lceopt(mus),lslack(mus));
    end

end

%%
jac = -jacobian(mus_len,qjac)';
jacf = matlabFunction(jac,'Vars',{q});
mus_lenf = matlabFunction(mus_len,'Vars',{q});
jacSize = size(jac);

if genEq == 1
    FQ = jac*mus_force';

    TEsc = invJtrans(q(5:8))*FQ(4:6);
    TEac = invJtrans(q(9:12))*FQ(7:9);
    TEgh = invJtrans(q(13:16))*FQ(10:12);
    TEel = FQ(13);

    TE_muscles = [zeros(13,1);TEsc;TEac;TEgh;TEel];
    matlabFunction(TE_muscles,'file','quaternion/TE_full_quat','vars',{t,q,act,fmax,lceopt,lslack});
    matlabFunction(mus_force,'file','quaternion/mus_forces_full_quat','vars',{t,q,act,fmax,lceopt,lslack});

end


% %%
% matlabFunction(mus_force,'file','mus_forces_simp_ulna','vars',{t,q,act,fmax,lceopt,lslack});
% matlabFunction(mus_len,'file','mus_lengths_simp_ulna','vars',{t,q});


end

function res = invJtrans(quat)
    q1 = quat(1);
    q2 = quat(2);
    q3 = quat(3);
    q4 = quat(4);
    res = [ q1/2,  q4/2, -q3/2;
            -q4/2,  q1/2,  q2/2;
            q3/2, -q2/2,  q1/2];
end
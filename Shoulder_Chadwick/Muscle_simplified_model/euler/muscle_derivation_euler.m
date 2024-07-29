function [jacf,mus_lenf,jacSize] = muscle_derivation_euler(model,genEq)
muscles = model.muscles;
% 
syms t real
q = sym('q', [14,1]);
act = sym('act', [length(muscles) 1]);
fmax = sym('fmax', [length(muscles) 1]);
lceopt = sym('lceopt', [length(muscles) 1]);
lslack = sym('lslack', [length(muscles) 1]);

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
                term = term * q(mdof);
            end
        end
        L = L + term;
    end
    mus_len(mus) = L;

    if genEq==1
        mus_force(mus) = muscle_force(act(mus),L,fmax(mus),lceopt(mus),lslack(mus));
    end

end
%%
jac = -jacobian(mus_len,q)';
jacf = matlabFunction(jac,'Vars',{q});
mus_lenf = matlabFunction(mus_len,'Vars',{q});
jacSize = size(jac);

if genEq==1
    FG = jac*mus_force';
    f_muscles = [zeros(10,1);FG(4:end-1)];
    %
    % matlabFunction(mus_force,'file','euler/mus_forces_simp_ulna','vars',{t,q,act,fmax,lceopt,lslack});
    % matlabFunction(mus_len,'file','euler/mus_lengths_simp_ulna','vars',{t,q});
    matlabFunction(f_muscles,'file','euler/FG_muscles_simp_ulna','vars',{t,q,act,fmax,lceopt,lslack});
end

end


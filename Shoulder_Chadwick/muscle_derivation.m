clear all
addpath("tmp\")
addpath("Functions\")
load("das3_full.mat")
muscles = model_full.muscles;

syms t real
q = sym('q', [10,1]);
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
                if mdof>10
                    term = term * 0;
                else
                    term = term * q(mdof);
                end
            end
        end
        L = L + term;
    end
    mus_len(mus) = L;
    mus_force(mus) = muscle_force(act(mus),L,fmax(mus),lceopt(mus),lslack(mus));

end
%%
mus_moment_arms = jacobian(mus_len,q);
f_muscles = [zeros(10,1);-mus_moment_arms'*mus_force'];
%%
matlabFunction(mus_force,'file','mus_forces_full_ulna','vars',{t,q,act,fmax,lceopt,lslack});
matlabFunction(mus_len,'file','mus_lengths_full_ulna','vars',{t,q});
matlabFunction(f_muscles,'file','f_muscles_full_ulna','vars',{t,q,act,fmax,lceopt,lslack});

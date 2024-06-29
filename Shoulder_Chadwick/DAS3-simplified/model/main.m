clear all
% addpath("build\")
% addpath("newpol\")
% model = das3_readosim('das3_simplified.osim','musclefile.mat','ghfile.mat');
% save('model_simplified','model')

addpath("build\tmp\")
addpath("newpol")
load("model_simplified.mat")
muscles = model.muscles;
%%
q = sym('q',[1,11]);
mus = randi([1,size(muscles,2)],1);
% mus = 90;
name = muscles{mus}.name
npolterms = muscles{mus}.lparam_count;
polcoeff = muscles{mus}.lcoefs;
expon = muscles{mus}.lparams;
musdof = muscles{mus}.dof_indeces-3;
nmusdof = muscles{mus}.dof_count;
L = 0; % Initialize the muscle length

for i = 1:npolterms
    % Add this term's contribution to the muscle length
    term = polcoeff(i);
    for j = 1:nmusdof
        mdof = musdof(j);
        for k = 1:expon(i, j)
            term = term * q(mdof); % This creates polcoeff[i] * product of all q_j to the power expon[i][j]
        end
    end
    L = L + term;
end
L
%%

load(sprintf('path_%s.mat',name))
%%

Lf = matlabFunction(L,'Vars',{q});
%%
n = randi([1,size(alljnts,1)],1);
% n = 50;
qs = alljnts(n,4:end);
% length_poly = Lf(qs(1),qs(2),qs(3),qs(4),qs(5),qs(6))
length_poly = Lf(qs)
length_os = alllengths(n)
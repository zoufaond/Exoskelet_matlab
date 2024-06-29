addpath('build\tmp\');
addpath('newpol\');
addpath('build\');
% osimfile = 'das3_simplified.osim';
% path = 'newpol';
% das3_polynomials(osimfile,path,'musclefile','ghfile')

model_simpl = das3_readosim('das3_simplified.osim','musclefile.mat');
save('das3_simplified','model_simpl')
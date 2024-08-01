addpath('polys_mod_simplified\')
addpath('build\');
osimfile = 'das3_mod_simplified.osim';
path = 'polys_mod_simplified';
das3_polynomials_quat(osimfile,path,'muscle_file_simplified_quat')

%%
model_simplified_quat = das3_readosim('das3_mod_simplified.osim','muscle_file_simplified_quat.mat');
save('das3_simplified_quat','model_simplified_quat')
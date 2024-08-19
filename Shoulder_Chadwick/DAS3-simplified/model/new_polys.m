addpath('polys_mod_full\')
addpath('build\');
osimfile = 'das3_mod_full.osim';
path = 'polys_mod_full';
das3_polynomials_quat(osimfile,path,'muscle_file_full_quat')

%%
model_full_quat = das3_readosim('das3_mod_full.osim','muscle_file_full_quat.mat');
save('das3_full_quat','model_full_quat')
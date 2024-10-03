addpath('polys_from_abd\')
addpath('build\');
%%
osimfile = 'das3_mod_full.osim';
path = 'polys_from_abd';
das3_polynomials_eul(osimfile,path,'muscle_file_eul')

%%

model_full_eul = das3_readosim('das3_mod_full.osim','muscle_file_eul.mat');
save('das3_eul_abd','model_full_eul')

%%
osimfile = 'das3_mod_full.osim';
path = 'polys_from_abd';
das3_polynomials_quat(osimfile,path,'muscle_file_quat')

%%

model_full_quat = das3_readosim('das3_mod_full.osim','muscle_file_quat.mat');
save('das3_quat_abd','model_full_quat')

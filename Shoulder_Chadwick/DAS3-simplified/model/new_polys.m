addpath('polys_mod_full_02\')
addpath('build\');
%%
osimfile = 'das3_mod_full.osim';
path = 'polys_mod_full_02';
das3_polynomials(osimfile,path,'muscle_file_full_eul_02')

%%

model_full_eul = das3_readosim('das3_mod_full.osim','muscle_file_full_eul_02.mat');
save('das3_full_eul_03','model_full_eul')

%%
osimfile = 'das3_mod_full.osim';
path = 'polys_mod_full_02';
das3_polynomials_quat(osimfile,path,'muscle_file_full_quat_02')

%%

model_full_quat = das3_readosim('das3_mod_full.osim','muscle_file_full_quat_02.mat');
save('das3_full_quat_03','model_full_quat')

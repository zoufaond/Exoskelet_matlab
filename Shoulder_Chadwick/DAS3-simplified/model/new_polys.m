addpath('newpol_modifiedModel\')
addpath('build\');
osimfile = 'das3_clav_scap_orig.osim';
path = 'newpol_modifiedModel';
das3_polynomials_quat(osimfile,path,'muscle_file_quat')

%%
model_full_quat = das3_readosim('das3_clav_scap_orig.osim','muscle_file_quat.mat');
save('das3_full_quat','model_full_quat')
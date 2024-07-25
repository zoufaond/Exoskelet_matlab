addpath('newpol_quat\')
addpath('build\');
osimfile = 'das3_simplified.osim';
path = 'newpol_quat';
das3_polynomials_quat(osimfile,path,'musclefile_quat')

model_simpl_quat = das3_readosim('das3_simplified.osim','musclefile_quat.mat');
save('das3_simplified_quaternion','model_simpl_quat')
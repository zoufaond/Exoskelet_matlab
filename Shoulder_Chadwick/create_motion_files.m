clearvars
motion_folder = 'Abduction';
osim_file = 'das3.osim';
mydir = 'Polyfiles';
musclepoly_file = 'musclepoly';

das3_polynomials(osim_file,mydir,musclepoly_file);
model = das3_readosim(osim_file,[mydir '/' musclepoly_file]);
save('OS_model','model')
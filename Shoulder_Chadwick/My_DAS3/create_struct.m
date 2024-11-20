clearvars
osim_file = 'das3.osim';
xml_file = strrep(osim_file,'.osim','.xml');
mydir = 'polyfiles';
musclepoly_file = 'musclepoly';

% das3_polynomials(osim_file,mydir,musclepoly_file);
model = das3_readosim(osim_file,[mydir '/' musclepoly_file]);
save('OS_model','model')
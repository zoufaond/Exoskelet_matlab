function data_to_mot(simdata,name)
% sim_data = load(simdata);
% angles = sim_data.angles;
angles = simdata;

fid = fopen(name,'w');
fprintf(fid, 'simulation\n');
fprintf(fid,'nRows=%i\n',size(angles,2));
fprintf(fid,'nColumns=%i\n',size(angles,1)+1);
fprintf(fid,'endheader\n');
fprintf(fid,'time  SC_y  SC_z  SC_x  AC_y  AC_z  AC_x  GH_y  GH_z  GH_yy  EL_x  PS_y\n');

for i=1:size(angles,2)
    for j=1:size(angles,1)
        if j==1
            fprintf(fid, '%4f', angles(j,i));
        elseif j == size(angles,1)
            fprintf(fid, '  %4f  0.000000\n', angles(j,i));
        else
            fprintf(fid, '  %4f', angles(j,i));
        end
    end
end
fclose(fid);
 % 
 % for i
 % if fid > 0
 %     fprintf(fid,'%d  %d  %d\n',myData');
 %     
 % end
end
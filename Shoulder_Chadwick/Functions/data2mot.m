function data2mot(out,name,type)


time = out.tout;

if type == 2
    angles = out.simdata_Q.signals.values(:,:)';
    SC_Q = angles(:,1:4);
    AC_Q = angles(:,5:8);
    GH_Q = angles(:,9:12);

    SC_RM = quat2rotm(SC_Q);
    AC_RM = quat2rotm(AC_Q);
    GH_RM = quat2rotm(GH_Q);

    SC_Eul = rotm2eul(SC_RM,'YZX');
    AC_Eul = rotm2eul(AC_RM,'YZX');
    GH_Eul = rotm2eul(GH_RM,'YZY');

    EL_x = angles(:,13);

    angles_OS = [SC_Eul,AC_Eul,GH_Eul,EL_x]*180/pi;

elseif type == 1
    angles = out.simdata_Eul.signals.values(:,:)';
    angles_OS = angles(:,1:10)*180/pi;
end

angles_OS = angles_OS';

fid = fopen(name,'w');
fprintf(fid, 'simulation\n');
fprintf(fid,'nRows=%i\n',size(angles_OS,2));
fprintf(fid,'nColumns=%i\n',size(angles_OS,1)+5);
fprintf(fid,'endheader\n');
fprintf(fid,'time  TH_x TH_y TH_z SC_y  SC_z  SC_x  AC_y  AC_z  AC_x  GH_y  GH_z  GH_yy  EL_x  PS_y\n');

for i=1:size(angles_OS,2)
    for j=1:size(angles_OS,1)
        if j==1
            fprintf(fid, '%4f  0.000000  0.000000  0.000000  %4f', time(i),angles_OS(j,i));
        elseif j == size(angles_OS,1)
            fprintf(fid, '  %4f  0.000000\n', angles_OS(j,i));
        else
            fprintf(fid, '  %4f', angles_OS(j,i));
        end
    end
end
fclose(fid);

end
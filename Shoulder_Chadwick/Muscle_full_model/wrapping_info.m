% clearvars
% 
% struct_os = readstruct('../DAS3-simplified/model/das3_clav_scap_orig_geommod_valmod.xml');
% 
% save('OS_struct.mat','struct_os')
clearvars
addpath ..\
OpenSimStruct = load('../OS_struct');


muscles = OpenSimStruct.struct_os.Model.ForceSet.objects.Schutte1993Muscle_Deprecated;

for imus = 1:2
    O_pos = muscles(imus).GeometryPath.PathPointSet.objects.PathPoint(1).location;
    I_pos = muscles(imus).GeometryPath.PathPointSet.objects.PathPoint(2).location;
    origin = muscles(imus).GeometryPath.PathPointSet.objects.PathPoint(1).socket_parent_frame;
    insertion = muscles(imus).GeometryPath.PathPointSet.objects.PathPoint(2).socket_parent_frame;
    isWrapped = isfield(muscles(imus).GeometryPath,'PathWrapSet');

    if strcmp(origin,'/bodyset/thorax') && isWrapped == 0 
        analytic(imus) = 1;
    else
        analytic(imus) = 0;
    end
end



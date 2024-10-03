addpath ../../../
addpath ../../../Motion/abduciton/
clearvars

modelEul = load("das3_eul.mat");
modelQuat = load("das3_quat.mat");

muscles = modelEul.model_full_eul.muscles;
folder = 'clavmod_IK';
dofs_names = {'SCy','SCz','SCx','ACy','ACz','ACx','GHy','GHz','GHyy','ELx','PSy'};
for i = 1:length(dofs_names)
    momarms.dof{i} = readtable([folder,'/',dofs_names{i},'_momarms.txt']);
end

lengths_tbl = readtable(['clavmod_IK/all_lengths.txt']);
motion = 'abd';
motion_file = load([motion,'_struct.mat']);

%%
for i = 1:length(muscles)
    imuscle = muscles{1,i};
    idofs = imuscle.dof_indeces;
    iname_table = imuscle.osim_name;
    iname = imuscle.name;
    for dof = 1:length(idofs)
        current_dof = idofs(dof)-3;
        imomarms = momarms.dof{1,current_dof};
        allmomarms(:,dof) = imomarms{:,iname_table};
    end
    alllengths = lengths_tbl{:,iname_table};
    alljnts = [zeros(151,3),motion_file.coords_struct.mot_euler_mod,ones(151,1)*1e-6];
    allmomarms(isnan(allmomarms)) = 1e-6;
    save(['path_',iname,'.mat'],'alljnts','allmomarms','alllengths')

    clear alljnts allmomarms alllengths
end



addpath Functions/

%choose two results struct to compare
folder_path = 'Motions/Scabduction_GL/';
motion_name = 'scabduction_GL.mat';

% folders = {[folder_path,'results_quat_QuatInit_100.mat']};
OS_struct = load([folder_path,motion_name]);
muscle_group = {'serr_ant'};
OS_model = [folder_path,'OS_model.mat'];

addpath([folder_path,'Poly_functions/'])

weights = {'150'}; %,'150','200','250'};
for iweight = 1:length(weights)
    folders = {[folder_path,'results_euler_QuatInit_',weights{iweight},'.mat'],[folder_path,'results_quat_QuatInit_',weights{iweight},'.mat']};
    plot_kinematics(folders,OS_struct);
    % plot_activations(folders,muscle_group,OS_model);
    plot_muscles(folders,muscle_group, OS_model);
end
% 
% 
% muscle_group = {'trap_scap','serr_ant','delt_scap','delt_clav','infra'};
% 
% EMG_struct = 'Experimental_data/EMG_struct.mat';
% plot_activations_EMG(folders, EMG_struct, OS_model);
% plot_activations(folders,muscle_group,OS_model)

function plot_activations_EMG(folders, EMG_struct, OS_model)

EMG_muscles = {'AnteriorDelt','IntermediateDelt','PosteriorDelt','Infrasp','MiddleTrap','UpperTrap','Serrupper', 'Serrlower'};
model_names = {'delt_scap11','delt_scap10','delt_scap_3','infra_1','trap_scap_7','trap_scap10','serr_ant_5','serr_ant_2'};

simulation = load(folders{2});
ending_val = 0;
time_simulation = simulation.data.tout(1:end-ending_val);

emg_data = load(EMG_struct);
model = load(OS_model);
muscles = model.model.muscles;
num_muscles = length(muscles);

for i = 1:num_muscles
    muscle_names{i} = muscles{i}.osim_name;
end


figure
tiledlayout(4,2);

current_emg_rsmpld = zeros(length(time_simulation),size(emg_data.data.(EMG_muscles{1}),2));
for i = 1:length(EMG_muscles)
    nexttile
    current_emg_full = emg_data.data.(EMG_muscles{i});
    current_emg = current_emg_full(1:end/2,:);
    time = linspace(0,1.5,size(current_emg,1));
    for ipar = 1:size(current_emg,2)
        current_emg_rsmpld(:,ipar) = spline(time,current_emg(:,ipar),time_simulation);
    end
    [S,M] = std(current_emg_rsmpld,0,2);
    upper_bound = M+S;
    lower_bound = M-S;
    plot(time_simulation,M,'k')
    hold on
    plot(time_simulation,upper_bound,'k',time_simulation,lower_bound,'k')
    hold on
    mus_index = find(strcmp(muscle_names,model_names{i}));
    plot(time_simulation,simulation.data.inputs(1:end-ending_val,mus_index),'r','LineWidth',1.5)
    hold on
    % patch([time' fliplr(time')], [lower_bound fliplr(upper_bound)], 'g')
    fill([time_simulation'; flip(time_simulation')],[lower_bound; flip(upper_bound)], 'b', 'edgecolor', 'none', 'facealpha', 0.1);
    title(EMG_muscles{i})
    axis([-inf inf 0 inf])
    xlabel('Time[s]')
    ylabel('Activation[-]')
    legend('Exp','','','Sim');
end

% fig = gcf;
% fig.Position(3) = fig.Position(3) + 250;
% Lgnd = legend('Exp','','','Sim');
% Lgnd.Position(1) = 0.9;
% Lgnd.Position(2) = 0.9;
% sgtitle('Activations', 'Interpreter', 'none')


end

function plot_activations(folders, muscle_group, OS_model)

model = load(OS_model);
muscles = model.model.muscles;
num_muscles = length(muscles);

for i = 1:num_muscles
    muscle_names{i} = muscles{i}.osim_name;
end
for imus_group = 1:length(muscle_group)
    current_mus_group = muscle_group{imus_group};
    mask = startsWith(muscle_names,current_mus_group);
    num_in_group = nnz(mask);
    plot_rows = ceil(num_in_group/3);
    current_names = muscle_names(mask);


    figure
    tiledlayout(plot_rows,3);
    for j = 1:num_in_group
        nexttile
        for i = 1:length(folders)
        current_struct = load(folders{i});
        activations = current_struct.data.inputs(:,mask);
        time = current_struct.data.tout;
        plot(time,activations(:,j),'LineWidth',1);hold on
        xlabel('time')
        ylabel('a[-]')
        end
        title(current_names(j), 'Interpreter', 'none')
    end
    fig = gcf;
    fig.Position(3) = fig.Position(3) + 250;
    Lgnd = legend('eul','quat');
    Lgnd.Position(1) = 0.01;
    Lgnd.Position(2) = 0.5;
    sgtitle('Activations', 'Interpreter', 'none')
end

end

function plot_kinematics(folders, OS_struct)
dofs_names = {'SCy','SCz','SCx','ACy','ACz','ACx','GHy','GHz','GHyy','ELx','PSy'};
num_coords = 10;
figure
tiledlayout(4,3);

for i = 1:num_coords
    nexttile
    for j = 1:length(folders)
        jmot = load(folders{j});
        time = jmot.data.tout';
        trajectory = jmot.data.trajectories;
        if size(trajectory,2) == 13
            trajectory_euler = quat2eul_motion(trajectory);
        else
            trajectory_euler = trajectory;
        end
        plot(time,trajectory_euler(:,i)*180/pi,'--','LineWidth',1.5)
        hold on
        title(dofs_names{i})
        % plot(time, trajectory_quat(:,i),'-.','LineWidth',1.5)
        % hold on
    end
    time = jmot.data.tout';
    OS_interp = spline(OS_struct.mot_struct.time,OS_struct.mot_struct.euler(:,i),time');
    %
    plot(time,OS_interp*180/pi,'g','LineWidth',1.5)

    xlabel('Time [s]')
    ylabel('Angle [deg]')
    % legend({'Euler','Quat','Experiment'})



end

fig = gcf;
fig.Position(3) = fig.Position(3) + 250;
Lgnd = legend('Euler','Quat','Experiment');
Lgnd.FontSize = 13;
Lgnd.Position(1) = 0.4;
Lgnd.Position(2) = 0.1;
end

function plot_muscles(folders,muscle_group, OS_model)

struct_euler = load(folders{1});
motion_euler = struct_euler.data.trajectories;
time = struct_euler.data.tout';
numdata = size(motion_euler,1);
struct_quat = load(folders{2});
motion_quat_orig = struct_quat.data.trajectories;
motion_quat_euler = quat2eul_motion(motion_quat_orig);

model = load(OS_model);
num_mus = length(model.model.muscles);
for i=1:num_mus
    all_muscle_names{i} = model.model.muscles{i}.osim_name;
end

mask = startsWith(all_muscle_names, muscle_group);
muscle_names = all_muscle_names(1,mask);
dofs_names = {'SCy','SCz','SCx','ACy','ACz','ACx','GHy','GHz','GHyy','ELx','PSy'};
alljoints = {'YZX','YZX','YZY'};
indexes = find(mask);

for imus = 1:length(muscle_names)
    current_mus = model.model.muscles{indexes(imus)};
    current_name = current_mus.osim_name
    dof_names = current_mus.dof_names;
    imomarms = zeros(11,101);
    JQuatInJEul = zeros(11,101);
    % size(motion_euler)
    % momarms_eul = examine_momarms(current_mus.Euler,dof_names,motion_euler);
    for iframe=1:numdata
        jac_eul = jacobiannoise_eul(0,motion_euler(iframe,:)');
        imomarms(:,iframe) = jac_eul(:,indexes(imus));
        jacinspat = JacInSpatnoise_quat(0,motion_quat_orig(iframe,:)');
        ijacinspat = jacinspat(:,indexes(imus));

         for j = 1:3
            JQuatInJEulCur = GeomJ(motion_euler(iframe,(j-1)*3+1:(j-1)*3+3),alljoints{j})*(ijacinspat((j-1)*3+1:(j-1)*3+3));
            JQuatInJEul((j-1)*3+1:(j-1)*3+3,iframe) = JQuatInJEulCur;
         end

         JQuatInJEul(10:11,iframe) = jacinspat(10:11);

    end

    [lengths_euler, dLdq_euler] = opensim_get_polyvalues(motion_euler, indexes(imus), current_mus.dof_indeces);
    [lengths_quat, dLdq_quat] = opensim_get_polyvalues(motion_quat_euler, indexes(imus), current_mus.dof_indeces);
    max_momarm = max([max(abs(dLdq_euler),[],'all'),max(abs(dLdq_quat),[],'all'),max(abs(JQuatInJEul),[],'all'),max(abs(imomarms),[],'all')]) + 0.01;
    figure
    
    for j=1:length(current_mus.dof_indeces)
        subplot(length(current_mus.dof_indeces),2,j*2-1)
        plot(time,dLdq_euler(:,j),time,imomarms(current_mus.dof_indeces(j)-3,:))
        axis([-inf inf -max_momarm max_momarm])
        subplot(length(current_mus.dof_indeces),2,j*2)
        plot(time,dLdq_quat(:,j),time,JQuatInJEul(current_mus.dof_indeces(j)-3,:))
        axis([-inf inf -max_momarm max_momarm])
    end
    % title(current_name)

end

end


function index = find_mus_index(model,mus_names)
    muscles = model.model_full_eul.muscles;
    nmus = model.model_full_eul.nMus;

    for i = 1:length(mus_names)
        for j = 1:nmus
            if strcmp(mus_names(i),muscles{j}.osim_name)
                index(i) = j;
                break
            end
        end
    end

end

function res = invJtrans(quat)
    q1 = quat(1);
    q2 = quat(2);
    q3 = quat(3);
    q4 = quat(4);
    res = [ q1/2,  q4/2, -q3/2;
            -q4/2,  q1/2,  q2/2;
            q3/2, -q2/2,  q1/2];
end

function res = GeomJ(phi,seq)
    s2 = sin(phi(2));
    s3 = sin(phi(3));
    c2 = cos(phi(2));
    c3 = cos(phi(3));
    if seq == 'YZX'
        res = [s2, c2*c3, -s3*c2;0,s3,c3;1,0,0];

    elseif seq == 'YZY'
        res = [s2*c3, c2, s2*s3; -s3, 0 ,c3; 0, 1, 0];
    end
end

function res = G(Q)
    Q0 = Q(1);
    Q1 = Q(2);
    Q2 = Q(3);
    Q3 = Q(4);
    res = [-Q1, Q0, Q3, -Q2;
            -Q2,-Q3, Q0, Q1;
            -Q3, Q2, -Q1, Q0];
end

function [lengths, minusdLdq] = opensim_get_polyvalues(angles, iMus, Dofs)
% This function calculates the length and -dL/dq of muscle "Mus" about dof set 
% "Dofs" at a given angle matrix "angles" of opensim model "Mod"
%
% "Angles" can be a vector (one hand position) or a matrix (one hand
% position per row)
% "Dofs" is a vector of indeces
%
% Adapted from opensim_get_momarm.m 
% Dimitra Blana, March 2012
%
% 28/3/2012: Use setValue (with active constraints) instead of set state
% 1/10/2014: Simplified how the muscles and dofs are accessed
%
% 11/25/19: Update by Derek Wolf: iMus is used to access the muscle and a
% check is used to determine if the name (with no underscores) matches the
% name in Mus

import org.opensim.modeling.*
osimfile = 'das3.osim';
Mod = Model(osimfile);

% this is needed to get the GH lines of action in the scapular frame
SimEn = SimbodyEngine();
SimEn.connectSimbodyEngineToModel(Mod);
groundbody = Mod.getBodySet().get('thorax');
scapulabody = Mod.getBodySet().get('scapula_r');
% initialize the system to get the initial state
state = Mod.initSystem;

Mod.equilibrateMuscles(state);

% If we only want one moment arm, put it in a cell
CoordSet = Mod.getCoordinateSet();
num_request_dofs = length(Dofs);

% get the muscle
% currentMuscle = Mod.getMuscles().get(Mus);
currentMuscle = Mod.getMuscles().get(iMus-1);
% if ~strcmp(fixname(char(currentMuscle.getName())),Mus)
%     error('Current muscle name is incorrect')
% end

% angles matrix: one position per row
angles = [zeros(size(angles,1),3) angles zeros(size(angles,1),1)];
[nrows,ncols] = size(angles);
nDofs = CoordSet.getSize;

if ncols~=nDofs
    if nrows~=nDofs
        errordlg('Angle matrix not the right size','Input Error');
        minusdLdq = [];
        return;
    else
        angles = angles';
    end
end

% initialise matrices for output
lengths = zeros(size(angles,1),1);
minusdLdq = zeros(size(angles,1),num_request_dofs);
quat_J = zeros(size(angles,1),num_request_dofs);

input_vector = Dofs';
num_elements = numel(input_vector);
remainder = mod(num_elements, 3);
if remainder ~= 0
    num_to_add = 3 - remainder;
    input_vector = [input_vector, zeros(1, num_to_add)];
end

dofs_block = reshape(input_vector, 3, []).';

joints = [];
ndofsblock = size(dofs_block,1);
for i=1:ndofsblock
    currentdofs = dofs_block(i,:);
    if currentdofs == [4,5,6]
        joints = [joints 2];
    elseif currentdofs == [7,8,9]
        joints = [joints 3];
    elseif currentdofs == [10,11,12]
        joints = [joints 4];
    end
end

if nargout>3
    GHfvecs = zeros(size(angles,1),3); 
end

for istep = 1:size(angles,1)
    if ~mod(istep,50)
        disp(['Muscle ',char(currentMuscle.getName()), ' - step ',...
            num2str(istep),' of ', num2str(size(angles,1))]);
    end
    
    % set dof values for this step
    for idof = 1:nDofs
        currentDof = CoordSet.get(idof-1);    
        currentDof.setValue(state,angles(istep,idof),1);
    end

    % get GH force line of action (if the muscle crosses GH)
    if nargout>3
        muspath = currentMuscle.getGeometryPath();
        fdarray = ArrayPointForceDirection();
        fvec2 = Vec3();
        muspath.getPointForceDirections(state,fdarray);
        scap_pt_index = -1;
        % find the "effective" muscle attachment on the scapula
        for ipt=1:fdarray.getSize-1
            body1 = char(fdarray.get(ipt-1).frame); %4.0 uses frames not bodies
            body2 = char(fdarray.get(ipt).frame);
            if strcmp(body1,'scapula_r')&&strcmp(body2,'humerus_r')
                scap_pt_index=ipt;
                break;
            elseif strcmp(body2,'scapula_r')&&strcmp(body1,'humerus_r')
                scap_pt_index=ipt-1;
                break;
            end
        end
        
        if scap_pt_index==-1
            % find the "effective" muscle attachment on the clavicle
            % instead
            for ipt=1:fdarray.getSize-1
                body1 = char(fdarray.get(ipt-1).frame);
                body2 = char(fdarray.get(ipt).frame);
                if strcmp(body1,'clavicle_r')&&strcmp(body2,'humerus_r')
                    scap_pt_index=ipt;
                    break;
                elseif strcmp(body2,'clavicle_r')&&strcmp(body1,'humerus_r')
                    scap_pt_index=ipt-1;
                    break;
                end
            end
        end

         if scap_pt_index==-1
            % find the "effective" muscle attachment on the thorax
            % instead
            for ipt=1:fdarray.getSize-1
                body1 = char(fdarray.get(ipt-1).frame);
                body2 = char(fdarray.get(ipt).frame);
                if strcmp(body1,'thorax')&&strcmp(body2,'humerus_r')
                    scap_pt_index=ipt;
                    break;
                elseif strcmp(body2,'thorax')&&strcmp(body1,'humerus_r')
                    scap_pt_index=ipt-1;
                    break;
                end
            end
         end

         if scap_pt_index==-1
            % find the "effective" muscle attachment on the ulna
            % instead
            for ipt=1:fdarray.getSize-1
                body1 = char(fdarray.get(ipt-1).frame);
                body2 = char(fdarray.get(ipt).frame);
                if strcmp(body1,'scapula_r')&&strcmp(body2,'ulna_r')
                    scap_pt_index=ipt;
                    break;
                elseif strcmp(body2,'scapula_r')&&strcmp(body1,'ulna_r')
                    scap_pt_index=ipt-1;
                    break;
                end
            end
         end
         
         if scap_pt_index==-1
            % find the "effective" muscle attachment on the radius
            % instead
            for ipt=1:fdarray.getSize-1
                body1 = char(fdarray.get(ipt-1).frame);
                body2 = char(fdarray.get(ipt).frame);
                if strcmp(body1,'scapula_r')&&strcmp(body2,'radius_r')
                    scap_pt_index=ipt;
                    break;
                elseif strcmp(body2,'scapula_r')&&strcmp(body1,'radius_r')
                    scap_pt_index=ipt-1;
                    break;
                end
            end
        end

        % calculate muscle force direction at that point (in global frame)
        scap_pt = fdarray.get(scap_pt_index);
        fvec = scap_pt.direction();
        
        % transform to the scapular coordinate frame
        SimEn.transform(state,groundbody,fvec,scapulabody,fvec2);
        GHfvecs(istep,1)=fvec2.get(0);
        GHfvecs(istep,2)=fvec2.get(1);
        GHfvecs(istep,3)=fvec2.get(2);
    end    
    
    % get length
    lengths(istep) = currentMuscle.getLength(state);
    % currentMuscle.computeMomentArm(state,'SCx')

    % get moment arms
    for idof = 1:num_request_dofs
        dofindex  = Dofs(idof)-1;
        currentDof = CoordSet.get(Dofs(idof)-1);
        minusdLdq(istep,idof) = currentMuscle.computeMomentArm(state,currentDof);
        
        % set dof to original value
        % currentDof.setValue(state,angles(istep,Dofs(idof)),1); 

        % input_vector = Dofs';
        % num_elements = numel(input_vector);
        % remainder = mod(num_elements, 3);
        % if remainder ~= 0
        %     num_to_add = 3 - remainder;
        %     input_vector = [input_vector, zeros(1, num_to_add)];
        % end
        % jnts2change_seq = {'YZX','YZX','YZX','YZY'};
        % dofs_block = reshape(input_vector, 3, []).';
        % 
        % joints = [];
        % ndofsblock = size(dofs_block,1);
        % for i=1:ndofsblock
        %     currentdofs = dofs_block(i,:);
        %     if currentdofs == [1,2,3]
        %         joints = [joints 2];
        %     elseif currentdofs == [4,5,6]
        %         joints = [joints 3];
        %     elseif currentdofs == [7,8,9]
        %         joints = [joints 4];
        %     end
        % end
        % 
        % 
        % % for i=1:numdata
        % for j=1:length(joints)
        %     if angles(istep,8) == 0
        %         angles(istep,8) = angles(istep,8) + 1e-3;
        %     end
        %     currentjnt = angles(istep,(joints(j)-1)*3+1:(joints(j)-1)*3+3);
        %     currentQ = quaternions(istep,(joints(j)-1)*4+1:(joints(j)-1)*4+4);
        %     currentmomarm = minusdLdq(istep,(j-1)*3+1:(j-1)*3+3);
        %     spat_J = Jeul2spat(currentjnt,currentmomarm,jnts2change_seq{joints(j)});
        %     quat_J(istep,(j-1)*3+1:(j-1)*3+3) = spatial2Jquat(currentQ,spat_J);
        % end
        % 
        % if ndofsblock > length(joints)
        %     % how many revolute joints
        %     nrevjnt = nnz(dofs_block(end,:));
        %     quat_J(:,end-nrevjnt+1:end) = minusdLdq(:,end-nrevjnt+1:end);
        % end
    end
end

end


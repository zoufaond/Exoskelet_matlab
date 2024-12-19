addpath Functions/

%choose two results struct to compare
folder_path = 'Motions/Steering/';
motion_name = 'Steering.mat';

% folders = {[folder_path,'results_quat_QuatInit_100.mat']};
OS_struct = load([folder_path,motion_name]);
muscle_group = {'delt_scap'};
OS_model = [folder_path,'OS_model.mat'];
addpath([folder_path,'Poly_functions/']);

% weights = {'300'};

% weights = {'100','150','200','250','300'};
for iweight = 1:length(weights)
    % folders = {[folder_path,'res_euler_Steering_',weights{iweight},'.mat'],[folder_path,'res_quat_Steering_',weights{iweight},'.mat']};
    folders = {[folder_path,'res_euler_Steering_250.mat'],[folder_path,'res_quat_Steering_250.mat']};

    % folders = {[folder_path,'res_euler_Scabduction_GL2_',weights{iweight},'.mat'],[folder_path,'res_quat_Scabduction_GL2_',weights{iweight},'.mat']};
    % folders = {[folder_path,'results_euler_QuatInit_',weights{iweight},'.mat'],[folder_path,'results_quat_QuatInit_',weights{iweight},'.mat']};
    plot_kinematics(folders,OS_struct);
    % plot_conoid_length(folders);
    % plot_SCx(folders);

    plot_activations(folders,muscle_group,OS_model);
    % plot_muscles(folders,muscle_group,OS_model);
end
% 
% 
% muscle_group = {'trap_scap','serr_ant','delt_scap','delt_clav','infra'};
% 
% EMG_struct = 'Experimental_data/EMG_struct.mat';
% plot_activations_EMG(folders, EMG_struct, OS_model);
% plot_activations(folders,muscle_group,OS_model)

function plot_SCx(folders)
    GH_names = {'GHy','GHz','GHyy'};
    num_coords = 10;
    euler_struct = load(folders{1});
    time = euler_struct.data.tout;
    traj_eul = euler_struct.data.trajectories;
    quat_struct = load(folders{2});
    traj_quat_orig = quat_struct.data.trajectories;
    traj_quat_eul = quat2eul_motion(traj_quat_orig);

    for i = 1:length(time)
        newSCx_eul(i,:) = min_conoid_length(traj_eul(i,:));
        length_eul(i,:) = conoid_length(traj_eul(i,4:6));
        newSCx_quat(i,:) = min_conoid_length(traj_quat_eul(i,:));
    end

    figure
    subplot(2,1,1)
    plot(time,rad2deg(newSCx_eul(:,3)),time,rad2deg(traj_eul(:,3)),'*')
    legend('Min length','simulation')
    subplot(2,1,2)
    plot(time,length_eul)
    title('eul')

    figure
    plot(time,rad2deg(newSCx_quat(:,3)),time,rad2deg(traj_quat_eul(:,3)),'*')
    legend('Min length','simulation')
    title('quat')
    

end

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
        if i == 1
            activations_eul_obj = sum(activations.^2,'all');
        end
        if i == 2
            activations_quat_obj = sum(activations.^2,'all');
        end
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
    text1 = annotation('textbox',[0.7,0.12,0.1,0.1],'string',{['Act Eul = ', num2str(activations_eul_obj),newline,'Act Quat = ', num2str(activations_quat_obj)]});
    sgtitle('Activations', 'Interpreter', 'none')
end

end

function plot_kinematics(folders, OS_struct)
dofs_names = {'SCy','SCz','SCx','Yscap in global','Zscap in global','Xscap in global','Yhum in global','Zhum in global','YYhum in global','ELx','PSy'};
GH_names = {'GHy','GHz','GHyy'};
num_coords = 10;
euler_struct = load(folders{1});
time = euler_struct.data.tout;
traj_eul = euler_struct.data.trajectories;
quat_struct = load(folders{2});
traj_quat_orig = quat_struct.data.trajectories;
traj_quat_eul = quat2eul_motion(traj_quat_orig);
traj_eul_obj = create_objective_traj_eul(traj_eul);
traj_quat_eul_obj = create_objective_traj_eul(traj_quat_eul);
OS_interp = interp1(OS_struct.mot_struct.time,OS_struct.mot_struct.euler,time',"spline");
OS_interp_obj = create_objective_traj_eul(OS_interp);
OS_interp_quat = eul2quat_motion(OS_interp);
OS_interp_quat_obj = create_objective_traj_quat(OS_interp_quat);
traj_quat_obj = create_objective_traj_quat(traj_quat_orig);

euler_act_obj = sum(euler_struct.data.inputs.^2,'all');
quat_act_obj = sum(quat_struct.data.inputs.^2,'all');

figure
idof = 1;
iplot = 1;
for i = 1:10
    subplot(5,3,iplot)
    plot(time,rad2deg(traj_eul_obj(:,idof)),'--','LineWidth',1.5)
    hold on
    plot(time,rad2deg(traj_quat_eul_obj(:,idof)),'-.','LineWidth',1.5)
    hold on
    title(dofs_names{idof})
    if i ~= 3
        plot(time,rad2deg(OS_interp_obj(:,idof)),'g','LineWidth',1.5)
    end
    axis([-inf inf -inf inf])  

    xlabel('Time [s]')
    ylabel('Angle [deg]')

    if idof == 9
        for GHdof = 1:3
        GHidof = GHdof+6;
        GHiplot = GHdof+9;
        subplot(5,3,GHiplot)
        plot(time,rad2deg(traj_eul(:,GHidof)),'--','LineWidth',1.5)
        hold on
        plot(time,rad2deg(traj_quat_eul(:,GHidof)),'-.','LineWidth',1.5)
        hold on
        title(GH_names{GHdof})
        plot(time,rad2deg(OS_interp(:,GHidof)),'g','LineWidth',1.5)
        if GHdof == 2
            yline(0,'r')
        end
        axis([-inf inf min(rad2deg(OS_interp(:,GHidof)))-10 max(rad2deg(OS_interp(:,GHidof)))+10])
    
        xlabel('Time [s]')
        ylabel('Angle [deg]')
        end
        iplot = iplot+3;
    end
    idof = idof+1;
    iplot = iplot+1;

end


fig = gcf;
fig.Position(3) = fig.Position(3) + 250;
Lgnd = legend('Euler','Quat','Experiment');
Lgnd.FontSize = 13;
Lgnd.Position(1) = 0.4;
Lgnd.Position(2) = 0.1;
text1 = annotation('textbox',[0.7,0.12,0.1,0.1],'string',{['Eul Act = ', num2str(euler_act_obj),newline,'Quat Act = ', num2str(quat_act_obj),]});

%%% QUAT KINEMATICS %%%

figure
idof = 1;
iplot = 1;
for i = 1:13
    subplot(5,4,i)
    plot(time,(traj_quat_obj(:,i)),'--','LineWidth',1.5)
    hold on
    plot(time,OS_interp_quat_obj(:,i))
    % title(dofs_names{idof})
    axis([-inf inf -inf inf])  

    xlabel('Time [s]')
    ylabel('Angle [deg]')

end

fig = gcf;
fig.Position(3) = fig.Position(3) + 250;
Lgnd = legend('Quat','Experiment');
Lgnd.FontSize = 13;
Lgnd.Position(1) = 0.4;
Lgnd.Position(2) = 0.1;

%%% END QUAT KINEMATICS %%%

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
alljoints = {'YZX','YZX','YZY'};
joints_names = {'SCy', 'SCz', 'SCx', 'ACy', 'ACz', 'ACx'};
indexes = find(mask);

for imus = 1:length(muscle_names)
    current_mus = model.model.muscles{indexes(imus)};
    current_name = current_mus.osim_name;
    dof_names = current_mus.dof_names;
    JQuatInJEul = zeros(numdata,11);
    dof_indeces = current_mus.dof_indeces-3;
    [lengths_eul,jacobian_eul] = momarms(current_mus.Euler, dof_indeces, motion_euler);
    motion_quat_WO_real = motion_quat_orig(:,[2:4,6:8,10:12,13]);
    [lengths_quat,jacobian_quat_current] = momarms(current_mus.Quaternion, dof_indeces, motion_quat_WO_real);
    jacobian_quat = zeros(numdata,11);
    jacobian_quat(:,dof_indeces) = jacobian_quat_current;

    for iframe=1:numdata
        for j = 1:3
            ind3 = ((j-1)*3+1:(j-1)*3+3);
            ind4 = ((j-1)*4+1:(j-1)*4+4);
            JQuatInSpat = invJtrans(motion_quat_orig(iframe,ind4)) * jacobian_quat(iframe,ind3)';
            JQuatInJEul(iframe,ind3) = GeomJ(motion_quat_euler(iframe,ind3),alljoints{j})*(JQuatInSpat);
        end
         JQuatInJEul(iframe,10) = jacobian_quat(10);
         JQuatInJEul(iframe,11) = jacobian_quat(11);
    end

    [lengths_OS_euler, dLdq_euler] = opensim_get_polyvalues(motion_euler, indexes(imus), current_mus.dof_indeces);
    [lengths_OS_quat, dLdq_quat] = opensim_get_polyvalues(motion_quat_euler, indexes(imus), current_mus.dof_indeces);
    max_momarm = max([max(abs(dLdq_euler),[],'all'),max(abs(dLdq_quat),[],'all'),max(abs(JQuatInJEul),[],'all'),max(abs(jacobian_eul),[],'all')]) + 0.01;
    figure
    
    for j=1:length(current_mus.dof_indeces)
        subplot(length(current_mus.dof_indeces),2,j*2-1)
        plot(time,dLdq_euler(:,j),time,-jacobian_eul(:,j))
        if j==1
            title(['Eul, Number of params:', string(current_mus.Euler.lparam_count)])
        end
        % plot(time,dLdq_euler(:,j),time,imomarms(current_mus.dof_indeces(j)-3,:))
        axis([-inf inf -max_momarm max_momarm])
        subplot(length(current_mus.dof_indeces),2,j*2)
        plot(time,dLdq_quat(:,j),time,-JQuatInJEul(:,dof_indeces(j)))
        if j==1
            title(['Quat, Number of params:', string(current_mus.Quaternion.lparam_count)'])
        end
        axis([-inf inf -max_momarm max_momarm])
    end
    % title(current_name)
    legend('osim','my')
    sgtitle(current_name)

    if strcmp(current_name,'serr_ant_6')

        figure
        tiledlayout(3,3);
        nexttile([1,3])
        plot(time,lengths_eul,'--','LineWidth',1.5);hold on
        plot(time,lengths_quat,'-.','LineWidth',1.5);hold on
        plot(time,lengths_OS_euler,'g','LineWidth',1.5);hold off
        title('Muscle length and moment arms approximation')
        ylabel('length [m]')
        xlabel('time [s]')
        axis([-inf inf 0 0.2])

        for j=1:length(current_mus.dof_indeces)
            nexttile
            plot(time,-jacobian_eul(:,j),'--','LineWidth',1.5);hold on
            plot(time,-JQuatInJEul(:,dof_indeces(j)),'-.','LineWidth',1.5);hold on
            plot(time,dLdq_euler(:,j),'g','LineWidth',1.5);
            axis([-inf inf -max_momarm max_momarm])
            ylabel('MA [m]')
            xlabel('time [s]')
            title(joints_names{j})
        end
        fig = gcf;
        fig.Position(3:4)=[700,350];
        Lgnd = legend('Euler approx','Quat approx','OpenSim');
        Lgnd.FontSize = 11;
        Lgnd.Position(1) = 0.7;
        Lgnd.Position(2) = 0.65;

        LRMS_eul = sqrt(sum((lengths_OS_euler-lengths_eul').^2,'all')/numdata)*1000
        LRMS_quat = sqrt(sum((lengths_OS_euler-lengths_quat').^2,'all')/numdata)*1000
        RRMS_eul = sqrt(sum((dLdq_euler+jacobian_eul(:,1:6)).^2,'all')/numdata)*1000
        RRMS_quat = sqrt(sum((dLdq_euler+JQuatInJEul(:,1:6)).^2,'all')/numdata)*1000

        exportgraphics(fig,'Serr_ant_12_approx.png','Resolution',600);
        
        

    end
end

% %%%%%%%%%%%%%%%%% plot with one trajectory %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% for imus = 1:3 %length(muscle_names)
%     current_mus = model.model.muscles{indexes(imus)};
%     current_name = current_mus.osim_name;
%     dof_names = current_mus.dof_names;
%     JQuatInJEul = zeros(numdata,11);
%     dof_indeces = current_mus.dof_indeces-3;
%     jacobian_eul = momarms(current_mus.Euler, dof_indeces, motion_quat_euler);
%     motion_quat_WO_real = motion_quat_orig(:,[2:4,6:8,10:12,13]);
%     jacobian_quat_current = momarms(current_mus.Quaternion, dof_indeces, motion_quat_WO_real);
%     jacobian_quat = zeros(numdata,11);
%     jacobian_quat(:,dof_indeces) = jacobian_quat_current;
% 
%     for iframe=1:numdata
%         for j = 1:3
%             ind3 = ((j-1)*3+1:(j-1)*3+3);
%             ind4 = ((j-1)*4+1:(j-1)*4+4);
%             JQuatInSpat = invJtrans(motion_quat_orig(iframe,ind4)) * jacobian_quat(iframe,ind3)';
%             JQuatInJEul(iframe,ind3) = GeomJ(motion_quat_euler(iframe,ind3),alljoints{j})*(JQuatInSpat);
%         end
%          JQuatInJEul(iframe,10) = jacobian_quat(10);
%          JQuatInJEul(iframe,11) = jacobian_quat(11);
%     end
% 
%     [lengths_osim, dLdq_osim] = opensim_get_polyvalues(motion_quat_euler, indexes(imus), current_mus.dof_indeces);
%     max_momarm = max([max(abs(dLdq_euler),[],'all'),max(abs(dLdq_quat),[],'all'),max(abs(JQuatInJEul),[],'all'),max(abs(jacobian_eul),[],'all')]) + 0.01;
% 
%     figure
%     for j=1:length(current_mus.dof_indeces)
%         subplot(floor(length(current_mus.dof_indeces)/3),3,j)
%         plot(time,-jacobian_eul(:,j))
%         hold on
%         plot(time,-JQuatInJEul(:,dof_indeces(j)))
%         hold on
%         plot(time,dLdq_osim(:,j))
% 
%         axis([-inf inf -max_momarm max_momarm])
%     end
%     legend('Euler','Quaternion','Osim')
    % sgtitle([current_name,' Eul, Quat nterms:',string(current_mus.Euler.lparam_count),string(current_mus.Quaternion.lparam_count)]) 

% end

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

function res = mulQuat(qa,qb)
    res = [ qa(1)*qb(1) - qa(2)*qb(2) - qa(3)*qb(3) - qa(4)*qb(4);
            qa(1)*qb(2) + qa(2)*qb(1) + qa(3)*qb(4) - qa(4)*qb(3);
            qa(1)*qb(3) - qa(2)*qb(4) + qa(3)*qb(1) + qa(4)*qb(2);
            qa(1)*qb(4) + qa(2)*qb(3) - qa(3)*qb(2) + qa(4)*qb(1)];
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

function rot_phix = R_x(phix)
    rot_phix = [1,0        , 0        ,0;
                0,cos(phix),-sin(phix),0;
                0,sin(phix), cos(phix),0;
                0,0        , 0        ,1];
end

function rot_phiy = R_y(phiy)
    rot_phiy = [cos(phiy),0,sin(phiy),0;
                0        ,1,0        ,0;
               -sin(phiy),0,cos(phiy),0;
                0        ,0,0        ,1];
end

function rot_phiz = R_z(phiz)
    rot_phiz = [cos(phiz),-sin(phiz),0,0;
                sin(phiz), cos(phiz),0,0;
                0           ,0      ,1,0;
                0           ,0      ,0,1];
end

function res = YZX_seq(angles)
    res = R_y(angles(1)) * R_z(angles(2)) * R_x(angles(3));
end

function res = YZY_seq(angles)
    res = R_y(angles(1)) * R_z(angles(2)) * R_y(angles(3));
end

function r = position(vec)
    r = [vec(1);vec(2);vec(3);1];
end

function trans = T_trans(vec)
    trans = [1,0,0,vec(1);
               0,1,0,vec(2);
               0,0,1,vec(3);
               0,0,0,1];
end

function res = Qrm(q)
    % rotation matrix from quaternion
    w = q(1);
    x = q(2);
    y = q(3);
    z = q(4);
    Rq =  [1-2*(y^2+z^2), 2*(x*y-z*w), 2*(x*z+y*w);
     2*(x*y+z*w), 1-2*(x^2+z^2), 2*(y*z-x*w);
     2*(x*z-y*w), 2*(y*z+x*w), 1-2*(x^2+y^2)];
    res = [Rq,zeros(3,1);
            zeros(1,3),1];
end

function res = create_objective_traj_eul(trajectory)
    res = zeros(size(trajectory));
    for istep = 1:size(trajectory,1)
        scapula_thorax = YZX_seq(trajectory(istep,1:3)) * YZX_seq(trajectory(istep,4:6));
        humerus_thorax = scapula_thorax * YZY_seq (trajectory(istep,7:9));
        res(istep,4:6) = rotm2eul(scapula_thorax(1:3,1:3),"YZX");
        res(istep,7:9) = rotm2eul(humerus_thorax(1:3,1:3),"YZY");
    end
    res(:,[1:3,10]) = trajectory(:,[1:3,10]);
end

function res = create_objective_traj_quat(trajectory)
    res = zeros(size(trajectory));
    for istep = 1:size(trajectory,1)
        scapula_thorax = mulQuat(trajectory(istep,1:4),trajectory(istep,5:8));
        humerus_thorax = mulQuat(scapula_thorax,trajectory(istep,(9:12)));
        res(istep,5:8) = scapula_thorax;
        res(istep,9:12) = humerus_thorax;
    end
    res(:,[1:4,13]) = trajectory(:,[1:4,13]);
end

% function res = hum_in_glob_quat(trajectory)
%     res = zeros(size(trajectory));
%     for istep = 1:size(trajectory,1)
%         scapula_thorax = mulQuat(trajectory(1:4),trajectory(5:8));
%         humerus_thorax = mulQuat(scapula_thorax,trajectory(9:12));
%         res(istep,4:6) = rotm2eul(scapula_thorax(1:3,1:3),"YZX");
%         res(istep,7:9) = rotm2eul(humerus_thorax(1:3,1:3),"YZY");
%     end
%     res(:,[1:3,10]) = trajectory(:,[1:3,10]);
% end


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
osimfile = 'das3_noserrant12.osim';
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

function [L,pmoment_arms] = momarms(musmodel, dof_indeces, angles)
% plot momentarm-angle data

% choose a subset of "angles" that contains only 100 points

% ...or use all angles
indeces = 1:size(angles,1);
sangles = angles(:,dof_indeces);

% calculate moment arms from polynomial
pmoment_arms = zeros(length(indeces),length(dof_indeces));
for iframe = 1:length(indeces)
    for i=1:musmodel.lparam_count

        % add this term's contribution to the muscle length 
        term = musmodel.lcoefs(i);

        for j=1:length(dof_indeces)
            for k=1:musmodel.lparams(i,j)
                term = term * sangles(iframe,j); % this creates lcoeff(i) * product of all angles to the power lparams(i,j) 
            end
        end

        % first derivatives of length with respect to all q's
        for  k=1:length(dof_indeces)
            % derivative with respect to q_k is zero unless exponent is 1 or higher and q is not zero
            if ((musmodel.lparams(i,k) > 0) && (sangles(iframe,k)))	
                dterm = musmodel.lparams(i,k)*term/sangles(iframe,k);
                pmoment_arms(iframe,k) = pmoment_arms(iframe,k) + dterm;
            end
        end
    end
end

L = zeros(1,length(indeces)); % Initialize the muscle length
for iframe = 1:length(indeces)
    Lterm = 0;
    for i=1:musmodel.lparam_count
        % Add this term's contribution to the muscle length
        term = musmodel.lcoefs(i);
        for j = 1:length(dof_indeces)
            for k = 1:musmodel.lparams(i, j)
                term = term * sangles(iframe,j);
            end
        end
        Lterm = Lterm + term;
    end

    L(iframe) = Lterm;
end

end

function res = min_conoid_length(traj)
    fun = @(x) conoid_length(x(2:4));
    % fun = @(x) sum((Rscap - R_scap_glob([traj(1),traj(2),x(1)],[x(2),x(3),x(4)])).^2,"all");
    nonlcon = @(x) mycon(x,traj);
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    lb = ones(1,4)*(-pi);
    ub = ones(1,4)*pi;
    x0 = traj(3:6);
    x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon);
    
    mot_eul_modified = traj;
    mot_eul_modified(3:6) = x;
    res = mot_eul_modified;
end

function res = conoid_length(AC)
    O = [0.1165,-0.0041,0.0143];
    I = T_trans([0.1575,0,0]) * YZX_seq(AC) * position([-0.0536, -0.0009, -0.0266]);
    res = sqrt((O(1) - I(1))^2 + (O(2) - I(2))^2 + (O(3) - I(3))^2);
end

function [c,ceq] = mycon(x,traj)
    Rscap = R_scap_glob(traj(1:3),traj(4:6));
    ceq = Rscap - R_scap_glob([traj(1),traj(2),x(1)],[x(2),x(3),x(4)]);
    c = [];
end

function res = R_scap_glob(jnt1,jnt2)
    res = YZX_seq(jnt1) * YZX_seq(jnt2);
end
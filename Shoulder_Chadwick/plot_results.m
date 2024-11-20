addpath Functions/

%choose two results struct to compare
folder_path = 'Motions/Abduction/';
motion_name = 'abduction.mat';
folders = {[folder_path,'results_euler_EulInit.mat'],[folder_path,'results_quat_InitQuat.mat']};
OS_struct = load([folder_path,motion_name]);
% plot_kinematics(folders,OS_struct)

muscle_group = {'trap_scap','serr_ant','delt_scap','delt_clav','infra'};
OS_model = [folder_path,'OS_model.mat'];
EMG_struct = 'Experimental_data/EMG_struct.mat';
% plot_activations_EMG(folders, EMG_struct, OS_model);
plot_activations(folders,muscle_group,OS_model)

function plot_activations_EMG(folders, EMG_struct, OS_model)

EMG_muscles = {'AnteriorDelt','IntermediateDelt','Infrasp','MiddleTrap','UpperTrap','PosteriorDelt','Serrupper', 'Serrlower'};
model_names = {'delt_clav_1','delt_scap10','infra_1','trap_scap_4','trap_scap10','delt_scap_3','serr_ant_6','serr_ant_1'};

emg_data = load(EMG_struct);
model = load(OS_model);
muscles = model.model.muscles;
num_muscles = length(muscles);

for i = 1:num_muscles
    muscle_names{i} = muscles{i}.osim_name;
end

for i = 1:length(model_names)
    find(strcmp(muscle_names,model_names{i}))
end
figure
tiledlayout(4,2);
for i = 1:length(EMG_muscles)
    nexttile
    current_emg = emg_data.data.(EMG_muscles{i});
    time = linspace(0,1.5,size(current_emg,1));
    [S,M] = std(current_emg,0,2);
    upper_bound = M+S;
    lower_bound = M-S;
    plot(time,M,'k')
    hold on
    plot(time,upper_bound,'k',time,lower_bound,'k')
    hold on
    % patch([time' fliplr(time')], [lower_bound fliplr(upper_bound)], 'g')
    fill([time'; flip(time')],[lower_bound; flip(upper_bound)], 'b', 'edgecolor', 'none', 'facealpha', 0.1);
    title(EMG_muscles{i})
end


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
        plot(time,trajectory_euler(:,i),'--','LineWidth',1.5)
        hold on
        % plot(time, trajectory_quat(:,i),'-.','LineWidth',1.5)
        % hold on
    end
    time = jmot.data.tout';
    OS_interp = spline(OS_struct.mot_struct.time,OS_struct.mot_struct.euler(:,i),time');
    %
    plot(time,OS_interp)



end

legend({'e','q','OS'})

end

% figure
% tiledlayout(4,3);
% for j = 1:num_coords
%     nexttile
%     OS_interp = spline(OS_struct.mot_struct.time,OS_struct.mot_struct.euler(:,i),time);
%     plot(time,trajectory_euler(:,i),'--', time, trajectory_quat(:,i),'-.',time,IK_interp,'LineWidth',1.5)
% end
% 


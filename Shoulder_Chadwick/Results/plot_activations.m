clearvars
addpath ../
model = load('das3_eul.mat');
muscles = model.model_full_eul.muscles;
num_muscles = length(muscles);
muscle_group = 'supra';
motion = 'abd'; %steering or abd

folders = {['eul_',motion],['quat_',motion]};

for i = 1:num_muscles
    muscle_names{i} = muscles{i}.osim_name;
end
mask = startsWith(muscle_names,muscle_group);
num_in_group = nnz(mask);
plot_rows = ceil(num_in_group/3);
current_names = muscle_names(mask);

figure
tiledlayout(plot_rows,3);
for j = 1:num_in_group
    nexttile
    for i = 1:length(folders)
    current_struct = load([folders{i},'/',folders{i},'.mat']);
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


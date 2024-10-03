addpath ../Functions/
addpath ../Motion
addpath ../Motion/abduciton/
clearvars

motion = 'abd'; %steering or abd

folders = {['eul_',motion,'_eulIC'],['quat_',motion,'_eulIC']};
dofs_names = {'SCy','SCz','SCx','ACy','ACz','ACx','GHy','GHz','GHyy','ELx','PSy'};
num_coords = 10;

struct1 = load([folders{1},'/',folders{1},'.mat']);
trajectory_euler = struct1.data.trajectories;
struct2 = load([folders{2},'/',folders{2},'.mat']);
trajectory_quat = quat2eul_motion(struct2.data.trajectories);
time = struct1.data.tout;
OS_motion = load([motion,'_struct.mat']);
trajectory_IK = OS_motion.coords_struct.mot_euler_mod;
time_IK = OS_motion.coords_struct.tout;



figure
tiledlayout(4,3);
for i = 1:num_coords
    nexttile
    IK_interp = spline(time_IK,trajectory_IK(:,i),time);
    plot(time,trajectory_euler(:,i)*180/pi, time, trajectory_quat(:,i)*180/pi,time,IK_interp*180/pi,'LineWidth',1)
    xlabel('time')
    ylabel('angle[°]')
    title(dofs_names{i})
end
fig = gcf;
fig.Position(3) = fig.Position(3) + 250;
Lgnd = legend('eul','quat','IK');
Lgnd.Position(1) = 0.01;
Lgnd.Position(2) = 0.5;
sgtitle('Trajectory', 'Interpreter', 'none')
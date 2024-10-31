clearvars
addpath ..\Muscle_full_model\euler\
addpath ..\Muscle_full_model\quaternion\
addpath ..\Muscle_full_model\motion_analysis\
addpath ..\
addpath ..\Motion\
addpath ..\Motion\abduciton\
addpath ..\Functions\
%%
modelEul = load("das3_eul.mat");
modelQuat = load("das3_quat.mat");
%%
% 'raw_IK/' or 'clavmod_IK/'
folder = 'clavmod_IK/';
muscle_group ='delt_scap';
motion = 'abd';

lengths_tbl = readtable(['motion_analysis/',folder,'all_lengths.txt']);
mask = startsWith(lengths_tbl.Properties.VariableNames, muscle_group);
muscle_names = lengths_tbl.Properties.VariableNames;
muscle_names = muscle_names(1,mask);
lengths_vals = lengths_tbl{:,mask};
motion_file = load([motion,'_struct.mat']);
motion_euler = motion_file.coords_struct.mot_euler_mod;
motion_quat = motion_file.coords_struct.mot_quaternion_mod;
numdata = length(motion_euler(:,1));
motion_euler = [motion_euler,zeros(numdata,1)];
motion_quat = [motion_quat,zeros(numdata,1)];
dofs_names = {'SCy','SCz','SCx','ACy','ACz','ACx','GHy','GHz','GHyy','ELx','PSy'};
alljoints = {'YZX','YZX','YZY'};

for i = 1:length(dofs_names)
    momarms.dof{i} = readtable(['motion_analysis/',folder,dofs_names{i},'_momarms.txt']);
end

time = motion_file.coords_struct.tout;
indexes = find_mus_index(modelEul,muscle_names);

% 
for imus = 10:10 %length(muscle_names)
    for i = 1:numdata
        LEul_current = mus_lengths_eul_abd(0,motion_euler(i,:)');
        iLEul(:,i) = LEul_current(indexes(imus));
        LQuat_current = mus_lengths_quat_abd(0,motion_quat(i,:)');
        iLQuat(:,i) = LQuat_current(indexes(imus));


        JEul_current = jacobian_eul_abd(0,motion_euler(i,:)');
        iJEul(:,i) = JEul_current(:,indexes(imus));

        JQuat_current = JacInSpat_quat_abd(0,motion_quat(i,:)');
        iJQuatInSpat = JQuat_current(:,indexes(imus));

        for j = 1:3
            JQuatInJEulCur = GeomJ(motion_euler(i,(j-1)*3+1:(j-1)*3+3),alljoints{j})*iJQuatInSpat((j-1)*3+1:(j-1)*3+3);
            JQuatInJEul((j-1)*3+1:(j-1)*3+3,i) = JQuatInJEulCur;
        end

        JQuatInJEul(10:11,i) = iJQuatInSpat(10:11);
    end

    RMS_length_quat = sqrt((sum((iLQuat.' - lengths_vals(:,imus)).^2))/numdata) * 1000;
    RMS_length_eul = sqrt((sum((iLEul.' - lengths_vals(:,imus)).^2))/numdata) * 1000;
    figure
    
    tiledlayout(2,3);
    nexttile([1,3])
    plot(time,iLEul,'g--',time,iLQuat,'r--',time,lengths_vals(:,imus),'b','LineWidth',1)
    title('Muscle length',['Quaternion RMS = ',num2str(RMS_length_quat,'%.3f'),'mm',', Euler RMS = ', num2str(RMS_length_eul,'%.3f'),'mm'], 'Interpreter', 'none')

    for iplot = 7:9 %length(dofs_names)
        momarmsvals = momarms.dof{iplot};
        momarmsvals = momarmsvals{:,mask};
        RMS_momarm_eul = sqrt((sum((iJEul(iplot,:).' - momarmsvals(:,imus)).^2))/numdata) * 1000;
        RMS_momarm_quat = sqrt((sum((JQuatInJEul(iplot,:).' - momarmsvals(:,imus)).^2))/numdata) * 1000;
        nexttile
        plot(time,iJEul(iplot,:)*1000,'g--',time,JQuatInJEul(iplot,:)*1000,'r--',time,momarmsvals(:,imus)*1000,'b','LineWidth',1.5) %
        title(string(dofs_names(iplot)) + ' moment arm',{['Quaternion RMS = ',num2str(RMS_momarm_quat,'%.3f'),'mm'],['Euler RMS = ', num2str(RMS_momarm_eul,'%.3f'),'mm']}, 'Interpreter', 'none')
        ylabel('Moment arm [mm]')
        axis([-inf inf -30 30])
    end
    % sgtitle(muscle_names(imus), 'Interpreter', 'none')
    fig = gcf;
    fig.Position(3) = fig.Position(3) + 250;
    Lgnd = legend('Euler','Quaternion','OpenSim');
    Lgnd.Position(1) = 0.8;
    Lgnd.Position(2) = 0.8;

end
exportgraphics(fig,'approximation_elevation.png','Resolution',600);

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

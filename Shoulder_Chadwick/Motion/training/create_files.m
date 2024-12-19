clearvars
addpath ..\..\Functions\
addpath ..\..\
motion_name = 'steering';
motion_file = readmatrix([motion_name,'.txt']);
motion_file(:,15) = scale_dof(motion_file(:,15),0.0);
[time_training, trajectory_training,~,~] = motion_polyfit(motion_file,'nplot_polyfits',40,40);
[~,~,time_simulation, trajectory_simulation] = motion_polyfit(motion_file,'nplot_polyfits',30,100);
numdata_training = length(trajectory_training);
osim_file = 'das3.osim';
OS_model = load("OS_model.mat");
data2mot(trajectory_simulation,time_simulation,[motion_name,'_simulation.txt'],'euler', 'struct');
data2mot(trajectory_simulation,time_simulation,[motion_name,'_simulation.mot'],'euler', 'struct');

%%


% for ival = 1:size(trajectory_simulation,1)
%     [~,Eeq1(ival),Eeq2(ival)] = is_close_to_thorax(trajectory_simulation(ival,:),OS_model,0,0);
% end
% plot(Eeq1)
% hold on
% plot(Eeq2)
%%
trajectory_training_full = [];
change_vals = linspace(-10,10,10)*pi/180;
for i=1:numdata_training
    
    if is_in_thorax_humerus(trajectory_training(i,:),OS_model) == 1 || is_in_thorax_scap(trajectory_training(i,:),OS_model) == 1
        continue
    end
    trajectory_training_min_con = min_conoid_length(trajectory_training(i,:));
    trajectory_training_full = [trajectory_training_full; trajectory_training_min_con]; 
    for ival = 1:length(change_vals)
        trajectory_rotclav_range = change_clavx(trajectory_training_min_con,trajectory_training_min_con(3) + change_vals(ival));
        trajectory_training_full = [trajectory_training_full; trajectory_rotclav_range]; 
    end
    % trajectory_training_full = [trajectory_training_full; trajectory_training(i,:)]; 

    % change_vals = linspace(-35,35,10)*pi/180;
    % for ival = 1:length(change_vals)
    %     trajectory_rotclav = change_clavx(trajectory_training(i,:),trajectory_training(i,3) + change_vals(ival));
    %     trajectory_training_full = [trajectory_training_full; trajectory_rotclav]; 

    num_noised = 30;
    noise_upper = [15,15,0,15,15,15,15,10,10,10,10];
    noise_lower =-[15,15,0,15,15,15,15,10,10,10,10];
    for inoise = 1:num_noised
        noise = noise_lower + (noise_upper-noise_lower).*rand(1,11);
        rnd_noise = noise*pi/180;
        trajectory_noised = trajectory_training(i,:) + rnd_noise;
        % trajectory_noised = trajectory_rotclav + rnd_noise;
        if is_in_thorax_humerus(trajectory_noised,OS_model) == 1 || is_in_thorax_scap(trajectory_noised,OS_model) == 1
            continue
        end
        %
        [is_scap_close,~] = is_close_to_thorax(trajectory_noised,OS_model,0.6,0.5);
        if  is_scap_close == 0
            continue
        end

        % trajectory_training_full = [trajectory_training_full; trajectory_noised]; 

        trajectory_rotclav = min_conoid_length(trajectory_noised);
        trajectory_training_full = [trajectory_training_full; trajectory_rotclav]; 
        
        for ival = 1:length(change_vals)
            trajectory_rotclav_range = change_clavx(trajectory_rotclav,trajectory_rotclav(3) + change_vals(ival));
            trajectory_training_full = [trajectory_training_full; trajectory_rotclav_range]; 
        end
    end

end

time_training_full = linspace(0,1,size(trajectory_training_full,1));
data2mot(trajectory_training_full,time_training_full,[motion_name,'_training.mot'],'euler', 'struct');
%%
% q = trajectory_training(7,:);
% res = is_in_thorax(q,OS_model);
% % close = is_close_to_thorax(q,OS_model,1)
% for idata = 1:size(trajectory,1)
%     [res(idata),val(idata)] = is_close_to_thorax(trajectory(idata,:),OS_model,0.9,0.25);
% end
% plot(val)
% hold on 
% plot(res)
    
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

function [c,ceq] = mycon(x,traj)
    Rscap = R_scap_glob(traj(1:3),traj(4:6));
    ceq = Rscap - R_scap_glob([traj(1),traj(2),x(1)],[x(2),x(3),x(4)]);
    c = [];
end

function res = conoid_length(AC)
    O = [0.1165,-0.0041,0.0143];
    I = T_trans([0.1575,0,0]) * YZX_seq(AC) * position([-0.0536, -0.0009, -0.0266]);
    res = sqrt((O(1) - I(1))^2 + (O(2) - I(2))^2 + (O(3) - I(3))^2);
end

function [res,Eeq1,Eeq2] = is_close_to_thorax(q,OS_model,mean_val,range_val)
    mus1 = OS_model.model.muscles{25};
    mus2 = OS_model.model.muscles{33};
    IP_mus1 = scapula_insertion_pos(q,mus1,OS_model);
    IP_mus2 = scapula_insertion_pos(q,mus2,OS_model);
    Epos = [0 -0.1486 0.0591];
    Edim = [0.147 0.2079 0.0944];
    Eeq1 = elips_eq(IP_mus1,Epos,Edim);
    Eeq2 = elips_eq(IP_mus2,Epos,Edim);
    val = sqrt(Eeq1^2+Eeq2^2);
    if abs(mean_val - val) < range_val
        res = 1;
    else
        res = 0;
    end
end

function res = is_in_thorax_scap(q,OS_model)

    Epos = [0 -0.1486 0.0591];
    Edim = [0.147 0.2079 0.0944];
    for i = 25:36
        current_mus = OS_model.model.muscles{i};
        IP_scap = scapula_insertion_pos(q,current_mus,OS_model);
        Eeq = elips_eq(IP_scap,Epos,Edim);
        res = 0;
        if Eeq <0
            res = 1;
            break
        end
    end
end

function res = is_in_thorax_humerus(q,OS_model)
    Epos_all = [0 -0.1486 0.0591;
            0.04 -0.19 0.05;
            0.05 -0.16 0.04];
    Edim_all = [0.147 0.2079 0.09;
            0.1 0.15 0.1;
            0.106 0.15 0.09];
    which_elips = [1,1,1,2,3,3];
    
    j = 1;
    for imus = 89:94
        current_mus = OS_model.model.muscles{imus};
        IP_hum = humerus_insertion_pos(q,current_mus,OS_model);
        Eeq = elips_eq(IP_hum,Epos_all(which_elips(j),:),Edim_all(which_elips(j),:));
        res = 0;
        if Eeq <0
            res = 1;
            break
        end
        j = j +1;
    end
end

function res = elips_eq(IP,Epos,Edim)
    res = ((IP(1)-Epos(1))/Edim(1))^2+((IP(2)-Epos(2))/Edim(2))^2+((IP(3)-Epos(3))/Edim(3))^2-1;
end

function pos = scapula_insertion_pos(q,muscle,OS_model)
    jnts = OS_model.model.joints;
    offset_thorax = jnts{1,2}.location;
    offset_clavicle = jnts{1,5}.location;
    insertion = muscle.origin_position;
    
    pos = T_trans(offset_thorax) * R_y(q(1)) * R_z(q(2)) * R_x(q(3)) * T_trans(offset_clavicle) * R_y(q(4)) * R_z(q(5)) * R_x(q(6)) * position(insertion);
end

function pos = humerus_insertion_pos(q,muscle,OS_model)
    jnts = OS_model.model.joints;
    offset_thorax = jnts{1,2}.location;
    offset_clavicle = jnts{1,5}.location;
    offset_scapula = jnts{1,8}.location;
    insertion = muscle.insertion_position;
    pos = T_trans(offset_thorax) * R_y(q(1)) * R_z(q(2)) * R_x(q(3)) * T_trans(offset_clavicle) * R_y(q(4)) * R_z(q(5)) * R_x(q(6)) *T_trans(offset_scapula) * R_y(q(7)) * R_z(q(8)) * R_y(q(9)) * position(insertion);
end

function T = T_trans(vec)
    T = [1,0,0,vec(1);
         0,1,0,vec(2);
         0,0,1,vec(3);
         0,0,0,1];
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

function r = position(vec)
    r = [vec(1);vec(2);vec(3);1];
end

function res = R_scap_glob(jnt1,jnt2)
    res = YZX_seq(jnt1) * YZX_seq(jnt2);
end

function res = YZX_seq(phi_vec)
    res = R_y(phi_vec(1)) * R_z(phi_vec(2)) * R_x(phi_vec(3));
end
clearvars
addpath ..\..\Functions\
addpath ..\..\
motion_file = readmatrix('scabduction.txt');
[time, trajectory] = motion_polyfit(motion_file,'no plot_polyfits',30);
numdata = length(time);
osim_file = 'das3.osim';

trajectory_training = [];
OS_model = load("OS_model.mat");

for i=1:numdata
    trajectory_training = [trajectory_training; trajectory(i,:)];

    change_vals = linspace(-10,20,6)*pi/180;
    for ival = 1:length(change_vals)
        trajectory_rotclav = change_clavx(trajectory(i,:),change_vals(ival));
        trajectory_training = [trajectory_training; trajectory_rotclav];

        num_noised = 50;
        noise_upper = [20,20,0,10,25,10,10,10,10,10,10];
        noise_lower =-[20,20,0,10,25,10,10,10,10,10,10];
        for inoise = 1:num_noised
            noise = noise_lower + (noise_upper-noise_lower).*rand(1,11);
            rnd_noise = noise*pi/180;
            trajectory_noised = trajectory_rotclav + rnd_noise;
            if is_in_thorax(trajectory_noised,OS_model) == 1
                continue
            end
            % 
            if is_close_to_thorax(trajectory_noised,OS_model,0.9,0.3) == 0
                continue
            end
            trajectory_training = [trajectory_training; trajectory_noised];
        end

    end
end

data2mot(trajectory_training,'scabduction_training.mot','euler', 'struct');
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
    

function [res,val] = is_close_to_thorax(q,OS_model,mean_val,range_val)
    mus1 = OS_model.model.muscles{25};
    mus2 = OS_model.model.muscles{33};
    Eeq1 = elips_eq(q,mus1,OS_model);
    Eeq2 = elips_eq(q,mus2,OS_model);
    val = sqrt(Eeq1^2+Eeq2^2);
    if abs(mean_val - val) < range_val
        res = 1;
    else
        res = 0;
    end
end

function res = is_in_thorax(q,OS_model)
    for i = 25:26
        current_mus = OS_model.model.muscles{i};
        Eeq = elips_eq(q,current_mus,OS_model);
        res = 0;
        if Eeq <0
            res = 1;
            break
        end
    end
end

function res = elips_eq(q,muscle,OS_model)
    Epos = [0 -0.1486 0.0591];
    Edim = [0.147 0.2079 0.0944];
    IP = scapula_insertion_pos(q,muscle,OS_model);
    res = ((IP(1)-Epos(1))/Edim(1))^2+((IP(2)-Epos(2))/Edim(2))^2+((IP(3)-Epos(3))/Edim(3))^2-1;
end

function pos = scapula_insertion_pos(q,muscle,OS_model)
    jnts = OS_model.model.joints;
    offset_thorax = jnts{1,2}.location;
    offset_clavicle = jnts{1,5}.location;
    insertion = muscle.origin_position;
    
    pos = T_trans(offset_thorax) * R_y(q(1)) * R_z(q(2)) * R_x(q(3)) * T_trans(offset_clavicle) * R_y(q(4)) * R_z(q(5)) * R_x(q(6)) * position(insertion);
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

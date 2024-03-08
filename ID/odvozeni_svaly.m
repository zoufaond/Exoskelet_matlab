clc;clear all
Glenohumeral = readtable('Glenohumeral_muscles.xlsx');
Scapulothoracic = readtable('Scapulothoracic_muscles.xlsx');
names_GH = table2array(Glenohumeral(:,1));
names_ST = table2array(Scapulothoracic(:,1));
origins_GH = table2array(Glenohumeral(:,5));
origin_coords_GH = coords2arr(table2array(Glenohumeral(:,6)));
insertions_GH = table2array(Glenohumeral(:,7));
insertion_coords_GH = coords2arr(table2array(Glenohumeral(:,8)));
origins_ST = table2array(Scapulothoracic(:,5));
origin_coords_ST = coords2arr(table2array(Scapulothoracic(:,6)));
insertions_ST = table2array(Scapulothoracic(:,7));
insertion_coords_ST = coords2arr(table2array(Scapulothoracic(:,8)));

N_GH = length(names_GH);
N_ST = length(names_ST);
q = sym('q',[1 9]);
dq = sym('dq', [1 9]);
akt_GH = sym('akt_GH',[1 N_GH]);
f0m_GH = sym('f0m_GH',[1 N_GH]);
l0m_GH = sym('l0m_GH', [1 N_GH]);
akt_ST = sym('akt_ST',[1 N_ST]);
f0m_ST = sym('f0m_ST',[1 N_ST]);
l0m_ST = sym('l0m_ST', [1 N_ST]);
t = sym('t');
for i = 1:N_GH
    muscle_lens_GH(i) = muscle_length(origins_GH{i},insertions_GH{i},origin_coords_GH(i,:),insertion_coords_GH(i,:),q);
    muscle_forces_GH(i) = muscle_force(muscle_lens_GH(i),f0m_GH(i),akt_GH(i),l0m_GH(i));
end
for i = 1:N_ST
    muscle_lens_ST(i) = muscle_length(origins_ST{i},insertions_ST{i},origin_coords_ST(i,:),insertion_coords_ST(i,:),q);
    muscle_forces_ST(i) = muscle_force(muscle_lens_ST(i),f0m_ST(i),akt_ST(i),l0m_ST(i));
end

fe = [zeros(9,1);jacobian(muscle_lens_ST,q)'*muscle_forces_ST'+jacobian(muscle_lens_GH,q)'*muscle_forces_GH'];
matlabFunction(fe,'file','f_muscles2','vars',{t,[q,dq],akt_GH,f0m_GH,l0m_GH,akt_ST,f0m_ST,l0m_ST});


function force = muscle_force(length, F_iso, akt, l0m)
    f_gauss = 0.25;
    force = (((length / l0m)^3) * exp(8 * length / l0m - 12.9) + (exp(-(length / l0m - 1)^2 / f_gauss)) * akt) * F_iso;
end

function length = muscle_length(origin, insertion, O_pos, I_pos, q)
    if strcmp(origin, 'Thorax') && strcmp(insertion, 'Clavicle')
        O = position(O_pos(1), O_pos(2), O_pos(3));
        T_t = [0.006325 0.00693 0.025465];
        TT_C = T_x(T_t(1))*T_y(T_t(2))*T_z(T_t(3));
        % I = TT_C*R_x(q(1)) * R_y(q(2)) * R_z(q(3)) * position(I_pos(1), I_pos(2), I_pos(3)); %123 rot
        I = TT_C*R_y(q(1)) * R_x(q(2)) * R_z(q(3)) * position(I_pos(1), I_pos(2), I_pos(3)); %213 rot
        
    elseif strcmp(origin, 'Thorax') && strcmp(insertion, 'Scapula')
        O = position(O_pos(1), O_pos(2), O_pos(3));
        T_t = [0.006325 0.00693 0.025465];
        T_c = [-0.01433 0.02007 0.135535];
        TT_C = T_x(T_t(1))*T_y(T_t(2))*T_z(T_t(3));
        RW_C = R_y(q(1)) * R_x(q(2)) * R_z(q(3));
        TC_S = T_x(T_c(1))*T_y(T_c(2))*T_z(T_c(3));
        ACJ_rot = R_x(-0.52) * R_y(0.52) * R_z(0);
        RC_S = R_z(q(4)) * R_y(q(5)) * R_x(q(6));
        ACJ_unrot = R_z(0) * R_y(-0.52) * R_x(0.52);
        I = TT_C * RW_C * TC_S * ACJ_rot * RC_S * ACJ_unrot  * position(I_pos(1), I_pos(2), I_pos(3)); %

    elseif strcmp(origin, 'Clavicle') && strcmp(insertion, 'Scapula')
        O =position(O_pos(1), O_pos(2), O_pos(3));
        TC_S = T_z(-10);
        RC_S = R_x(q(4)) * R_y(q(5)) * R_z(q(6));
        I =TC_S * ACJ_rot * RC_S * ACJ_unrot * position(I_pos(1), I_pos(2), I_pos(3));

    elseif strcmp(origin, 'Clavicle') && strcmp(insertion, 'Humerus')
        T_c = [-0.01433 0.02007 0.135535];
        T_s = [-0.00955 -0.034 0.009];
        O =position(O_pos(1), O_pos(2), O_pos(3));
        TC_S = T_x(T_c(1))*T_y(T_c(2))*T_z(T_c(3));
        ACJ_rot = R_x(-0.52) * R_y(0.52) * R_z(0);
        RC_S = R_z(q(4)) * R_y(q(5)) * R_x(q(6));
        ACJ_unrot = R_z(0) * R_y(-0.52) * R_x(0.52);
        TS_H = T_x(T_s(1))*T_y(T_s(2))*T_z(T_s(3));
        RS_H = R_x(q(7)) * R_y(q(8)) * R_z(q(9));
        I =TC_S * ACJ_rot * RC_S * ACJ_unrot * TS_H * RS_H * position(I_pos(1), I_pos(2), I_pos(3));

    elseif strcmp(origin, 'Scapula') && strcmp(insertion, 'Humerus')
        T_s = [-0.00955 -0.034 0.009];
        O =position(O_pos(1), O_pos(2), O_pos(3));
        TS_H = T_x(T_s(1))*T_y(T_s(2))*T_z(T_s(3));
        RS_H = R_x(q(7)) * R_y(q(8)) * R_z(q(9));
        I =TS_H*RS_H*position(I_pos(1), I_pos(2), I_pos(3));

    elseif strcmp(origin, 'Thorax') && strcmp(insertion, 'Humerus')
        T_c = [-0.01433 0.02007 0.135535];
        T_s = [-0.00955 -0.034 0.009];
        T_t = [0.006325 0.00693 0.025465];
        TT_C = T_x(T_t(1))*T_y(T_t(2))*T_z(T_t(3));
        O =position(O_pos(1), O_pos(2), O_pos(3));
        RT_C = TT_C*R_y(q(1)) * R_x(q(2)) * R_z(q(3));
        TC_S = T_x(T_c(1))*T_y(T_c(2))*T_z(T_c(3));
        ACJ_rot = R_x(-0.52) * R_y(0.52) * R_z(0);
        RC_S = R_z(q(4)) * R_y(q(5)) * R_x(q(6));
        ACJ_unrot = R_z(0) * R_y(-0.52) * R_x(0.52);
        TS_H = T_x(T_s(1))*T_y(T_s(2))*T_z(T_s(3));
        RS_H = R_x(q(7)) * R_y(q(8)) * R_z(q(9));
        I =RT_C * TC_S * ACJ_rot * RC_S * ACJ_unrot * TS_H * RS_H * position(I_pos(1), I_pos(2), I_pos(3));

    end

    length = sqrt((O(1) - I(1))^2 + (O(2) - I(2))^2 + (O(3) - I(3))^2);
end

function trans_x = T_x(x)
    trans_x = [1,0,0,x;
               0,1,0,0;
               0,0,1,0;
               0,0,0,1];
end

function trans_y = T_y(y)
    trans_y = [1,0,0,0;
               0,1,0,y;
               0,0,1,0;
               0,0,0,1];
end

function trans_z = T_z(z)
    trans_z = [1,0,0,0;
               0,1,0,0;
               0,0,1,z;
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

function r = position(x,y,z)
    r = [x;y;z;1];
end

function t = coords2arr(coords)
    t = zeros(length(coords),3);
    for i = 1:length(coords)
        t(i,:) = str2num(coords{i});
    end
end
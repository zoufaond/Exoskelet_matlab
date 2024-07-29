function [EulTraj,QuatTraj] = create_abduction_traj(initCond,n,t_end,scale)

% a sine function for abduction trajectory
Abduction = (sin(linspace(0,pi,n)-pi/2)'+1)*scale;

    for i=1:n
        % this creates trajectory w.r.t. the scapula frame, scapula is
        % fixed
        qeuli = [initCond(1:7);initCond(8)+Abduction(i);initCond(9:10)];
        qeul(:,i) = qeuli;
        
        % this creates trajectory w.r.t. the thorax frame
        % Rhum = Humerus_to_Thorax_matrix(qn); 
        % [x(i),y(i),z(i)] = rotxyz(Rhum);
    end

% analytic differentation
% dqs = [zeros(7,n);scale*(cos(linspace(0,pi,n)-pi/2)*pi/t_end);zeros(2,n)];
% ddqs = [zeros(7,n);scale*(-sin(linspace(0,pi,n)-pi/2)*(pi/t_end).^2);zeros(2,n)];

% numerical differentiation
stepSize = t_end/n;
dqeul = diff(qeul,1,2)/stepSize;
ddqeul = diff(dqeul,1,2)/stepSize;

% create struct with euler angles
EulTraj.q = qeul;
EulTraj.dq = [zeros(10,1),dqeul];
EulTraj.ddq = [zeros(10,1),ddqeul,zeros(10,1)];


% joints
joints = {'YZX','YZX','YZY'};
allquat = zeros(13,n);

% transform euler angles to quaternion
for i=1:n %n
    for idof = 1:length(joints)
        euli = EulTraj.q((idof-1)*3+1:(idof-1)*3+3,i);
        quati = eul2quat(euli',joints{idof});
        allquat((idof-1)*4+1:(idof-1)*4+4,i) = quati;
    end
end

% last row is revolute joint, so it is the same
allquat(end,:) = qeul(end,:);
dallquat = diff(allquat,1,2)/stepSize;

% preallocate
allomega = zeros(10,n-1);

% map quaternion to angular velocity (omega)
for i = 1:n-1
    for idof = 1:length(joints)
        omegai = 2*G(allquat((idof-1)*4+1:(idof-1)*4+4,i))*dallquat((idof-1)*4+1:(idof-1)*4+4,i);
        allomega((idof-1)*3+1:(idof-1)*3+3,i) = omegai;
    end
end

allomega(end,:) = dqeul(end,:);
dallomega = diff(allomega,1,2)/stepSize;

QuatTraj.quat = allquat;
QuatTraj.dquat = [zeros(13,1),dallquat];
QuatTraj.omega = [zeros(10,1),allomega];
QuatTraj.domega = [zeros(10,1),dallomega,zeros(10,1)];
% traj = [x;y;z];
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
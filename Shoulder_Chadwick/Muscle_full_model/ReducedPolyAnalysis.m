function ReducedPolyAnalysis()

% this calculates the RMS for lengths and jacobians for quaternion and
% euler angles

addpath euler\
addpath quaternion\
addpath ..\
addpath ..\DAS3-simplified\model\polys_mod_full\
modelEul = load("das3_eul.mat");
modelQuat = load("das3_quat.mat");

% generate jacobians and muscle length (as anonymous functions that depends on q1..qn coordinates)
genEq = 0;
[JEul,LEul,EulSize] = muscle_derivation_euler(modelEul.model_full_eul,genEq);
[JQuat,LQuat,QuatSize] = muscle_derivation_quat(modelQuat.model_full_quat,genEq);

%%
% this is used only for muscle names
muscles = modelEul.model_full_eul.muscles;

% iterate over muscles
for imus = 1:length(muscles)
musfilename = ['..\DAS3-simplified\model\polys_mod_full\path_',muscles{imus}.name,'.mat'];
ma = load(musfilename);

mus = muscles{imus};
numdata = size(ma.alljnts,1);
alljoints = {'YZX','YZX','YZX','YZY'}; % list of joints that are transformed from Eul to Quat and vice versa
dofs = mus.dof_indeces;
ndofs = length(mus.dof_indeces);

% here we create the jacobians out of the momarms
alljac = zeros(size(ma.alljnts));

for i = 1:numdata
    alljac(i,dofs(1):dofs(end)) = ma.allmomarms(i,:);
end


for i = 1:numdata
    % get muscle lengths and jacobian of Eul polynomials
    JEulVal = JEul(ma.alljnts(i,:)');
    JEulPol(i,:) = JEulVal(:,imus)';

    LEulVal = LEul(ma.alljnts(i,:)');
    LEulPol(i,1) = LEulVal(imus);

    % get muscle lengths from Quat polynomials
    LQuatVal = LQuat(ma.alljntsQ(i,:)');
    LQuatPol(i,1) = LQuatVal(imus);

    % get jacobians from Quat polynomials
    JQuatVal = JQuat(ma.alljntsQ(i,:)');
    JQuatCur = JQuatVal(:,imus);

    % map from quaternion jacobians to euler jacobians (each joint separately)
    for j = 1:4
        JQuatInSpatCur = invJtrans(ma.alljntsQ(i,(j-1)*4+1:(j-1)*4+4))*JQuatCur((j-1)*3+1:(j-1)*3+3);
        JQuatInJEulCur = GeomJ(ma.alljnts(i,(j-1)*3+1:(j-1)*3+3),alljoints{j})*JQuatInSpatCur;
        JQuatInJEul(i,(j-1)*3+1:(j-1)*3+3) = JQuatInJEulCur';
    end

    % revolute joints are the same
    JQuatInJEul(i,13:14) = JQuatCur(13:14);
end

% compare calculated length from polynomials with the lengths exported from
% OpenSim and calculate RMS
resLEul = ma.alllengths-LEulPol;
resLQuat = ma.alllengths-LQuatPol;
RMSLEul = (sqrt(sum(resLEul.^2)/length(resLEul)));
RMSLQuat = (sqrt(sum(resLQuat.^2)/length(resLQuat)));

% compare Eul jacobians and Quat jacobians (mapped to the Eul space) and
% calculate RMS
resJEul = alljac-JEulPol;
resJQuat = alljac-JQuatInJEul;
RMSJEul = sqrt(sum(resJEul.^2,'all')/length(resJEul))*1000; 
RMSJQuat = sqrt(sum(resJQuat.^2,'all')/length(resJQuat))*1000; 

%save RMS for bar plots
RMSlen(imus,:) = [RMSLEul,RMSLQuat];
RMSJac(imus,:) = [RMSJEul,RMSJQuat];
xbar{imus} = mus.name;

clear LEulPol JEulPol LQuatPol JQuatInJEul resLEul resLQuat RMSLEul RMSLQuat resJEul resJQuat RMSJEul RMSJQuat
 
end

Sxbar = categorical(string(xbar));

figure
bar(Sxbar,RMSlen)
title('RMSlength')
legend('Eul','Quat')
set(gca, 'YScale', 'log')

figure
bar(Sxbar,RMSJac)
title('Jac')
legend('Eul','Quat')
set(gca, 'YScale', 'log')

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
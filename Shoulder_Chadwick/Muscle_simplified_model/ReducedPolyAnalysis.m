function ReducedPolyAnalysis()

addpath euler\
addpath quaternion\
addpath ..\
addpath ..\DAS3-simplified\model\newpol_quat\
modelEul = load("das3_simplified.mat");
modelQuat = load("das3_simplified_quaternion.mat");

genEq = 0;
[JEul,LEul,EulSize] = muscle_derivation_euler(modelEul.model_simpl,genEq);
[JQuat,LQuat,QuatSize] = muscle_derivation_quat(modelQuat.model_simpl_quat,genEq);

%%
muscles = modelEul.model_simpl.muscles;

for imus = 1:length(muscles)
musfilename = ['..\DAS3-simplified\model\newpol_quat\path_',muscles{imus}.name,'.mat'];
ma = load(musfilename);

mus = muscles{imus};
numdata = size(ma.alljnts,1);
alljoints = {'YZX','YZX','YZX','YZY'}; % joints to change
dofs = mus.dof_indeces;
ndofs = length(mus.dof_indeces);

alljac = zeros(size(ma.alljnts));

% create jacobian of momarms from OpenSim
for i = 1:numdata
    alljac(i,dofs(1):dofs(end)) = ma.allmomarms(i,:);
end



for i = 1:numdata
    JEulVal = JEul(ma.alljnts(i,:)');
    LEulVal = LEul(ma.alljnts(i,:)');
    LEulPol(i,1) = LEulVal(imus);
    JEulPol(i,:) = JEulVal(:,imus)';

    JQuatVal = JQuat(ma.alljntsQ(i,:)');
    LQuatVal = LQuat(ma.alljntsQ(i,:)');
    LQuatPol(i,1) = LQuatVal(imus);
    JQuatCur = JQuatVal(:,imus);

    for j = 1:4
        JQuatInSpatCur = invJtrans(ma.alljntsQ(i,(j-1)*4+1:(j-1)*4+4))*JQuatCur((j-1)*3+1:(j-1)*3+3);
        JQuatInJEulCur = GeomJ(ma.alljnts(i,(j-1)*3+1:(j-1)*3+3),alljoints{j})*JQuatInSpatCur;
        JQuatInJEul(i,(j-1)*3+1:(j-1)*3+3) = JQuatInJEulCur';
    end
    JQuatInJEul(i,13:14) = JQuatCur(13:14);
end

resLEul = ma.alllengths-LEulPol;
resLQuat = ma.alllengths-LQuatPol;
RMSLEul = (sqrt(sum(resLEul.^2)/length(resLEul)));
RMSLQuat = (sqrt(sum(resLQuat.^2)/length(resLQuat)));

resJEul = alljac-JEulPol;
resJQuat = alljac-JQuatInJEul;
RMSJEul = (sqrt(sum(resJEul.^2,'all')))/(numdata*ndofs); 
RMSJQuat = (sqrt(sum(resJQuat.^2,'all')))/(numdata*ndofs); 

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
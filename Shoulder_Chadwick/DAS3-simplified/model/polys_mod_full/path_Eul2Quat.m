clearvars
addpath ../
clc;
osimfile = '../das3_clav_scap_orig.osim';
model = das3_readosim(osimfile);
mydir = 'newpol_modifiedModel';
muscles = model.muscles;

for currentmuscle = 1:length(muscles)
imus = currentmuscle;
% 
musfilename = ['../newpol_modifiedModel/path_',muscles{imus}.name,'.mat'];
ma = load(musfilename);

mus = muscles{imus};

num_data = size(ma.alljnts,1);
fprintf(1,'Muscle name:      %s\n',mus.name);

input_vector = mus.dof_indeces';
num_elements = numel(input_vector);
remainder = mod(num_elements, 3);
if remainder ~= 0
    num_to_add = 3 - remainder;
    input_vector = [input_vector, zeros(1, num_to_add)];
end

dofs = reshape(input_vector, 3, []).';

alljoints = {'YZX','YZX','YZX','YZY'}; % joints to change
nalljoints = length(alljoints);
numdata = size(ma.alljnts,1);
%
alljntsQ = [];
alljntsQA = [];

% change angles seq to quaternions
for i=1:nalljoints
    quat = eul2quat(ma.alljnts(:,(i-1)*3+1:(i-1)*3+3),alljoints{i});
    alljntsQ = [alljntsQ quat]; % full quaternions
    alljntsQA = [alljntsQA quat(:,2:4)]; %only 2:4 term of quaternion (input to A matrix)
end

alljntsQ = [alljntsQ ma.alljnts(:,end-1:end)]; % add revolute joints
alljntsQA = [alljntsQA ma.alljnts(:,end-1:end)]; % add revolute joints

joints = [];
ndofsblock = size(dofs,1);
for i=1:ndofsblock
    currentdofs = dofs(i,:);
    if currentdofs == [4,5,6]
        joints = [joints 2];
    elseif currentdofs == [7,8,9]
        joints = [joints 3];
    elseif currentdofs == [10,11,12]
        joints = [joints 4];
    end
end

quat_J = zeros(size(ma.allmomarms));

for i=1:numdata
    for j=1:length(joints)
        currentjnt = ma.alljnts(i,(joints(j)-1)*3+1:(joints(j)-1)*3+3);
        currentQ = alljntsQ(i,(joints(j)-1)*4+1:(joints(j)-1)*4+4);
        currentmomarm = ma.allmomarms(i,(j-1)*3+1:(j-1)*3+3);
        spat_J = Jeul2spat(currentjnt,currentmomarm,alljoints{joints(j)});
        quat_J(i,(j-1)*3+1:(j-1)*3+3) = spatial2Jquat(currentQ,spat_J);
    end
end

if ndofsblock > length(joints)
    % how many revolute joints
    nrevjnt = nnz(dofs(end,:));
    quat_J(:,end-nrevjnt+1:end) = ma.allmomarms(:,end-nrevjnt+1:end);
else
    % do nothing
end

save(['path_',mus.name,'.mat'],'alljntsQ','alljntsQA','quat_J','-append')
end


function res = Jeul2spat(phi,EulJ,seq)
    s2 = sin(phi(2));
    s3 = sin(phi(3));
    c2 = cos(phi(2));
    c3 = cos(phi(3));
    if seq == 'YZX'
        res = [s2, c2*c3, -s3*c2;0,s3,c3;1,0,0]\EulJ';

    elseif seq == 'YZY'
        res = [s2*c3, c2, s2*s3; -s3, 0 ,c3; 0, 1, 0]\EulJ';
    end

    res = res';
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

function res = spatial2Jquat(quat,jac)
    res = Jt(quat) * jac';
    res = res';
end

function res = T(quat)
    a = quat(1);
    b = quat(2);
    c = quat(3);
    d = quat(4);
    res = [-b/a, -c/a, -d/a;eye(3)];
end

function res = Jt(quat)
    res = 2*G(quat)*T(quat);
    res = res';
end
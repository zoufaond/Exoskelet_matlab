function [JacInSpatf,mus_lenf,jacSize] = muscle_derivation_quat(model,genEq)
muscles = model.muscles;

syms t real
q = sym('q', [18,1]);
qpol = [q(2:4);q(6:8);q(10:12);q(14:16);q(17);q(18)];
act = sym('act', [length(muscles) 1]);
fmax = sym('fmax', [length(muscles) 1]);
lceopt = sym('lceopt', [length(muscles) 1]);
lslack = sym('lslack', [length(muscles) 1]);
% 
for mus=1:length(muscles)

    isWrapped = model.muscles{mus}.isWrapped;

    if isWrapped == 0
        origin = model.muscles{mus}.origin_frame;
        insertion = model.muscles{mus}.insertion_frame;
        O_pos = model.muscles{mus}.origin_position;
        I_pos = model.muscles{mus}.insertion_position;
        L = analytic_length(origin, insertion, O_pos, I_pos, q(5:end), model);
        Jac = -jacobian(L,q)';

        for dof = 1:4
        JacInSpat((dof-1)*3+1:(dof-1)*3+3,mus) = 1/2 * G(q((dof-1)*4+1:(dof-1)*4+4)) * Jac((dof-1)*4+1:(dof-1)*4+4);
        end
        JacInSpat(13:14,mus) = Jac(17:18);
        disp('analytic mus')

    else
        name = muscles{mus}.name;
        npolterms = muscles{mus}.lparam_count;
        polcoeff = muscles{mus}.lcoefs;
        expon = muscles{mus}.lparams;
        musdof = muscles{mus}.dof_indeces;
        nmusdof = muscles{mus}.dof_count;
        L = 0; % Initialize the muscle length
    
        for i = 1:npolterms
            % Add this term's contribution to the muscle length
            term = polcoeff(i);
            for j = 1:nmusdof
                mdof = musdof(j);
                for k = 1:expon(i, j)
                    term = term * qpol(mdof);
                end
            end
            L = L + term;
        end

        isWrappedVec(mus) = 1;
        Jac = -jacobian(L,qpol)';
        for dof = 1:4
            JacInSpat((dof-1)*3+1:(dof-1)*3+3,mus) = invJtrans(q((dof-1)*4+1:(dof-1)*4+4)) * Jac((dof-1)*3+1:(dof-1)*3+3);
        end
        JacInSpat(13:14,mus) = Jac(13:14);
        disp('poly mus')

    end

    mus_len(mus) = L;

    if genEq == 1
        mus_force(mus) = muscle_force(act(mus),L,fmax(mus),lceopt(mus),lslack(mus));
    end

end
disp('Jac done')
%%
JacInSpatf = matlabFunction(JacInSpat,'Vars',{q});
mus_lenf = matlabFunction(mus_len,'Vars',{q});
jacSize = size(JacInSpat);
disp('matlabfunctions done')

if genEq == 1
    FQ = JacInSpat*mus_force';
    disp('FQ done')

    TE_muscles = [zeros(13,1);FQ(4:end-1)];

    % matlabFunction(TE_muscles,'file','quaternion/TE_full_quat','vars',{t,q,act,fmax,lceopt,lslack});
    matlabFunction(mus_force,'file','quaternion/mus_forces_full_quat','vars',{t,q,act,fmax,lceopt,lslack});

end


% %%
% matlabFunction(mus_force,'file','mus_forces_simp_ulna','vars',{t,q,act,fmax,lceopt,lslack});
% matlabFunction(mus_len,'file','mus_lengths_simp_ulna','vars',{t,q});


end

function muscle_length = analytic_length(origin, insertion, O_pos, I_pos, q, model)
offset_thorax = model.joints{1,2}.location;
offset_clavicle = model.joints{1,5}.location;

    if strcmp(origin, 'thorax') && strcmp(insertion, 'clavicle_r')
        O = position(O_pos(1), O_pos(2), O_pos(3));
        I = T_trans(offset_thorax) * Qrm(q(1:4)) * position(I_pos(1), I_pos(2), I_pos(3));
        
    elseif strcmp(origin, 'thorax') && strcmp(insertion, 'scapula_r')
        O = position(O_pos(1), O_pos(2), O_pos(3));
        RW_C = T_trans(offset_thorax) * Qrm(q(1:4));
        TC_S = T_trans(offset_clavicle);
        RC_S = Qrm(q(5:8));
        I = RW_C * TC_S * RC_S * position(I_pos(1), I_pos(2), I_pos(3));
    end

    muscle_length = sqrt((O(1) - I(1))^2 + (O(2) - I(2))^2 + (O(3) - I(3))^2);
end

function res = Qrm(q)
    % rotation matrix from quaternion
    w = q(1);
    x = q(2);
    y = q(3);
    z = q(4);
    Rq =  [1-2*(y^2+z^2), 2*(x*y-z*w), 2*(x*z+y*w);
     2*(x*y+z*w), 1-2*(x^2+z^2), 2*(y*z-x*w);
     2*(x*z-y*w), 2*(y*z+x*w), 1-2*(x^2+y^2)];
    res = [Rq,sym(zeros(3,1));
            sym(zeros(1,3)),1];
end

function trans_x = T_trans(vec)
    trans_x = [1,0,0,vec(1);
               0,1,0,vec(2);
               0,0,1,vec(3);
               0,0,0,1];
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

function r = position(x,y,z)
    r = [x;y;z;1];
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
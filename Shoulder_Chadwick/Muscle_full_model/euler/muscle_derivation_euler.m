function [jacf,mus_lenf,jacSize] = muscle_derivation_euler(model,genEq)
muscles = model.muscles;

% 
syms t real
q = sym('q', [14,1]);
act = sym('act', [length(muscles) 1]);
fmax = sym('fmax', [length(muscles) 1]);
lceopt = sym('lceopt', [length(muscles) 1]);
lslack = sym('lslack', [length(muscles) 1]);

for mus=1:length(muscles)
    isWrapped = model.muscles{mus}.isWrapped;

    if isWrapped == 0
        origin = model.muscles{mus}.origin_frame;
        insertion = model.muscles{mus}.insertion_frame;
        O_pos = model.muscles{mus}.origin_position;
        I_pos = model.muscles{mus}.insertion_position;
        
        L = analytic_length(origin, insertion, O_pos, I_pos, q(4:end), model);
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
                    term = term * q(mdof);
                end
            end
            L = L + term;
        end
        disp('poly mus')

    end


    mus_len(mus) = L;

    if genEq==1
        mus_force(mus) = muscle_force(act(mus),L,fmax(mus),lceopt(mus),lslack(mus));
    end

end
%%
jac = -jacobian(mus_len,q)';
jacf = matlabFunction(jac,'Vars',{q});
mus_lenf = matlabFunction(mus_len,'Vars',{q});
jacSize = size(jac);

if genEq==1
    FG = jac*mus_force';
    f_muscles = [zeros(10,1);FG(4:end-1)];
    %
    matlabFunction(mus_force,'file','euler/mus_forces_full_eul','vars',{t,q,act,fmax,lceopt,lslack});
    matlabFunction(f_muscles,'file','euler/TE_full_eul','vars',{t,q,act,fmax,lceopt,lslack});
end

end

function muscle_length = analytic_length(origin, insertion, O_pos, I_pos, q, model)
offset_thorax = model.joints{1,2}.location;
offset_clavicle = model.joints{1,5}.location;

    if strcmp(origin, 'thorax') && strcmp(insertion, 'clavicle_r')
        O = position(O_pos(1), O_pos(2), O_pos(3));
        I = T_trans(offset_thorax) * R_y(q(1)) * R_z(q(2)) * R_x(q(3)) * position(I_pos(1), I_pos(2), I_pos(3));
        
    elseif strcmp(origin, 'thorax') && strcmp(insertion, 'scapula_r')
        O = position(O_pos(1), O_pos(2), O_pos(3));
        RW_C = T_trans(offset_thorax) * R_y(q(1)) * R_z(q(2)) * R_x(q(3));
        TC_S = T_trans(offset_clavicle);
        RC_S = R_y(q(4)) * R_z(q(5)) * R_x(q(6));
        I = RW_C * TC_S * RC_S * position(I_pos(1), I_pos(2), I_pos(3));
    end

    muscle_length = sqrt((O(1) - I(1))^2 + (O(2) - I(2))^2 + (O(3) - I(3))^2);
end

function trans_x = T_trans(vec)
    trans_x = [1,0,0,vec(1);
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

function r = position(x,y,z)
    r = [x;y;z;1];
end



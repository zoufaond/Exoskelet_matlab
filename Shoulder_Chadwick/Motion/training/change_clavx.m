function res = change_clavx(q, change_value)

clavicle_xrot_change = [0,0,change_value];
mot_eul_modified = zeros(size(q));

Rscap = R_scap_glob(q(1:3),q(4:6));
new_xrot = q(1:3) + clavicle_xrot_change;
fun = @(x) sum((Rscap - R_scap_glob(new_xrot,[x(1),x(2),x(3)])).^2,"all");
x0 = [q(4:6)];
A = [];
b = [];
Aeq = [];
beq = [];
x = fmincon(fun,x0,A,b,Aeq,beq);

mot_eul_modified(1:3) = new_xrot;
mot_eul_modified(4:6) = x;

mot_eul_modified(:,7:end) = q(:,7:end);

res = mot_eul_modified;

end


function res = R_scap_glob(jnt1,jnt2)
    res = YZX_seq(jnt1) * YZX_seq(jnt2);
end

function res = YZX_seq(phi_vec)
    res = R_y(phi_vec(1)) * R_z(phi_vec(2)) * R_x(phi_vec(3));
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
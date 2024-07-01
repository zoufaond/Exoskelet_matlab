function [traj,q,dq,ddq] = create_abduction_traj(initCond,n,t_end,scale)
Abduction = (sin(linspace(0,pi,n)-pi/2)'+1)*scale;
    for i=1:n
        qn = [initCond(1:7);initCond(8)+Abduction(i);initCond(9:10)];
        q(:,i) = qn;
        Rhum = Humerus_to_Thorax_matrix(qn);
        [x(i),y(i),z(i)] = rotxyz(Rhum);
    end
h = t_end/n;
dq = diff(q,1,2)/h;
dq = [dq,dq(:,end)];
ddq = diff(q,2,2)/(h^2);
ddq = [ddq(:,1),ddq,ddq(:,end)];
traj = [x;y;z];
end
function [traj,q,dqs,ddqs] = create_abduction_traj(initCond,n,t_end,scale)
Abduction = (sin(linspace(0,pi,n)-pi/2)'+1)*scale;
    for i=1:n
        qn = [initCond(1:7);initCond(8)+Abduction(i);initCond(9:10)];
        q(:,i) = qn;
        Rhum = Humerus_to_Thorax_matrix(qn);
        [x(i),y(i),z(i)] = rotxyz(Rhum);
    end
% h = t_end/n;
dqs = [zeros(7,n);scale*(cos(linspace(0,pi,n)-pi/2)*pi/t_end);zeros(2,n)];
ddqs = [zeros(7,n);scale*(-sin(linspace(0,pi,n)-pi/2)*(pi/t_end).^2);zeros(2,n)];
% dq = diff(q,1,2)/h;
% dq = [zeros(10,1),dq];
% ddq = diff(q,2,2)/(h^2);
% ddq = [ddq,ddq(:,end),ddq(:,end)];

traj = [x;y;z];
end
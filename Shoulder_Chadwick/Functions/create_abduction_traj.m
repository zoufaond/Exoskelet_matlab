function traj = create_abduction_traj(initCond,n)
scale = 0.8;
Abduction = (sin(linspace(0,pi,n)-pi/2)'+1)*scale;
    for i=1:n
        Rhum = Humerus_to_Thorax_matrix([initCond(1:7);initCond(8)+Abduction(i);initCond(9)]);
        [x(i),y(i),z(i)] = rotxyz(Rhum);
    end
traj = [x;y;z];
end
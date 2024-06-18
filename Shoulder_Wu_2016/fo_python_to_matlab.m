syms t x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16 x17 x18 m1 m2 m3 c k g IMxx IMyy IMzz IMxy IMyz IMzx IUxx IUyy IUzz IUxy IUyz IUzx ILxx ILyy ILzz ILxy ILyz ILzx rigid0Px rigid0Py rigid0Pz rigid2Cx rigid2Cy rigid2Cz rigid1Px rigid1Py rigid1Pz rigid1Cx rigid1Cy rigid1Cz x_cont1 y_cont1 z_cont1 x_cont2 y_cont2 z_cont2 k_contact mx my mz ax ay az eps real
q = [x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16 x17 x18];
m = [m1,m2,m3];
IU = [IUxx IUyy IUzz IUxy IUyz IUzx];
IM = [IMxx IMyy IMzz IMxy IMyz IMzx];
IL = [ILxx ILyy ILzz ILxy ILyz ILzx];
rigid0P = [rigid0Px, rigid0Py, rigid0Pz];
rigid1C = [rigid1Cx, rigid1Cy, rigid1Cz];
rigid1P = [rigid1Px, rigid1Py, rigid1Pz];
rigid2C = [rigid2Cx, rigid2Cy, rigid2Cz];
rigid2P = [rigid2Px, rigid2Py, rigid2Pz];
rigid3C = [rigid3Cx, rigid3Cy, rigid3Cz];
cont_params = [mx,my,mz,ax,ay,az,x_cont1,y_cont1,z_cont1,x_cont2,y_cont2,z_cont2,k_contact, eps];



matlabFunction(fo,'file','fo','vars',{t,q,IU,IM,IL,m,c,k,g,rigid0P,rigid1C,rigid1P,rigid2C,rigid2P,rigid3C,cont_params});
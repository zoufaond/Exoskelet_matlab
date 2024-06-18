syms t x1 x2 x3 x4 x5 x6 x7 x8 x9 c k g rigid3Cx rigid3Cy rigid3Cz rigid2Px rigid2Py rigid2Pz rigid2Cx rigid2Cy rigid2Cz rigid1Px rigid1Py rigid1Pz rigid1Cx rigid1Cy rigid1Cz m1 m2 m3 IUxx IUyy IUzz IUxy IUyz IUzx IMxx IMyy IMzz IMxy IMyz IMzx ILxx ILyy ILzz ILxy ILyz ILzx real
q = [x1,x2,x3,x4,x5,x6,x7,x8,x9];
m = [m1,m2,m3];
IU = [IUxx IUyy IUzz IUxy IUyz IUzx];
IM = [IMxx IMyy IMzz IMxy IMyz IMzx];
IL = [ILxx ILyy ILzz ILxy ILyz ILzx];
rigid1C = [rigid1Cx, rigid1Cy, rigid1Cz];
rigid1P = [rigid1Px, rigid1Py, rigid1Pz];
rigid2C = [rigid2Cx, rigid2Cy, rigid2Cz];
rigid2P = [rigid2Px, rigid2Py, rigid2Pz];
rigid3C = [rigid3Cx, rigid3Cy, rigid3Cz];


matlabFunction(mm,'file','mm','vars',{t,q,IU,IM,IL,m,c,k,g,rigid1C,rigid1P,rigid2C,rigid2P,rigid3C});
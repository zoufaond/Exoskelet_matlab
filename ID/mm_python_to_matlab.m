syms t x1 x2 x3 x4 x5 x6 x7 x8 x9 c_m s_m h_m cixx ciyy cizz cixy ciyz cizx sixx siyy sizz sixy siyz sizx hixx hiyy hizz hixy hiyz hizx T_cx T_cy T_cz T_sx T_sy T_sz ccomx ccomy ccomz scomx scomy scomz hcomx hcomy hcomz
q = [x1,x2,x3,x4,x5,x6,x7,x8,x9];
m = [c_m,s_m,h_m];
cI = [cixx, ciyy, cizz, cixy, ciyz, cizx];
sI = [sixx, siyy, sizz, sixy, siyz, sizx];
hI = [hixx, hiyy, hizz, hixy, hiyz, hizx];
ccom = [ccomx, ccomy, ccomz];
scom = [scomx, scomy, scomz];
hcom = [hcomx, hcomy, hcomz];
T_c = [T_cx,T_cy,T_cz];
T_s = [T_sx,T_sy,T_sz];

% 'mm' z https://colab.research.google.com/drive/13ojY-aWEMXYoFUJRt8fr1WP0Z2mWlCNu#scrollTo=gMq-hb7mfPfT

matlabFunction(mm,'file','mm_wu_unrot','vars',{t,q,m,cI,sI,hI,ccom,scom,hcom,T_c,T_s});
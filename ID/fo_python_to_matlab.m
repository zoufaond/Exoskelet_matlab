syms t c k g x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16 x17 x18 c_m s_m h_m cixx ciyy cizz cixy ciyz cizx sixx siyy sizz sixy siyz sizx hixx hiyy hizz hixy hiyz hizx T_tx T_ty T_tz T_cx T_cy T_cz T_sx T_sy T_sz ccomx ccomy ccomz scomx scomy scomz hcomx hcomy hcomz x_cont1 y_cont1 z_cont1 x_cont2 y_cont2 z_cont2 k_contact mx my mz ax ay az eps
q = [x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16 x17 x18];
m = [c_m,s_m,h_m];
Ic = [cixx, ciyy, cizz, cixy, ciyz, cizx];
Is = [sixx, siyy, sizz, sixy, siyz, sizx];
Ih = [hixx, hiyy, hizz, hixy, hiyz, hizx];
ccom = [ccomx, ccomy, ccomz];
scom = [scomx, scomy, scomz];
hcom = [hcomx, hcomy, hcomz];
m_el = [mx, my, mz];
a_el = [ax,ay,az];
cont1 = [x_cont1,y_cont1,z_cont1];
cont2 = [x_cont2,y_cont2,z_cont2];
cont_params = [k_contact, eps];
T_t = [T_tx, T_ty, T_tz];
T_c = [T_cx,T_cy,T_cz];
T_s = [T_sx,T_sy,T_sz];
p = [c,k,g];


matlabFunction(fo,'file','fo_wu_unrot','vars',{t,q,m,Ic,Is,Ih,ccom,scom,hcom,T_t,T_c,T_s,p,m_el,a_el,cont1,cont2,cont_params});
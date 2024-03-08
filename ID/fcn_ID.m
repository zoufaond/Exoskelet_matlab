function torq = fcn_ID(t,q,dq,ddq,m,cI,sI,hI,ccom,scom,hcom,T_t,T_c,T_s,p,m_el,a_el,cont1,cont2,cont_params,akt_GH, Glenohumeral_F0M, Glenohumeral_l0,akt_ST, Scapulothoracic_F0M, Scapulothoracic_l0)
ccI = [cI(1:3),cI(6),cI(4),cI(5)];
ssI = [sI(1:3),sI(6),sI(4),sI(5)];
hhI = [hI(1:3),hI(6),hI(4),hI(5)];

torq = zeros(18,1);

fosim   = fo_wu(0,[[q(1),q(2),q(3),q(4),q(5),q(6),q(7),q(8),q(9)],[dq(1),dq(2),dq(3),dq(4),dq(5),dq(6),dq(7),dq(8),dq(9)]],m,cI,sI,hI,ccom,scom,hcom,T_t,T_c,T_s,p,m_el,a_el,cont1,cont2,cont_params);
fesim = f_muscles2(0,[q(1),q(2),q(3),q(4),q(5),q(6),q(7),q(8),q(9)],akt_GH, Glenohumeral_F0M, Glenohumeral_l0,akt_ST, Scapulothoracic_F0M, Scapulothoracic_l0);

mmsim = mm_wu(0,[q(1),q(2),q(3),q(4),q(5),q(6),q(7),q(8),q(9)],m,ccI,ssI,hhI,ccom,scom,hcom,T_c,T_s);
ddqq = [zeros(1,9),ddq(1),ddq(2),ddq(3),ddq(4),ddq(5),ddq(6),ddq(7),ddq(8),ddq(9)]';
torq = mmsim*ddqq-fosim-fesim;
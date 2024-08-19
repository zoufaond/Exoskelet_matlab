function res = quat2eul_motion(angles)

SC_Q = angles(:,1:4);
AC_Q = angles(:,5:8);
GH_Q = angles(:,9:12);

SC_RM = quat2rotm(SC_Q);
AC_RM = quat2rotm(AC_Q);
GH_RM = quat2rotm(GH_Q);

SC_Eul = rotm2eul(SC_RM,'YZX');
AC_Eul = rotm2eul(AC_RM,'YZX');
GH_Eul = rotm2eul(GH_RM,'YZY');

EL_x = angles(:,13);

res = [SC_Eul,AC_Eul,GH_Eul,EL_x];
end
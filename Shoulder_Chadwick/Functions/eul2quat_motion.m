function res = eul2quat_motion(angles)

SC_EUL = angles(:,1:3);
AC_EUL = angles(:,4:6);
GH_EUL = angles(:,7:9);

SC_RM = quat2rotm(SC_EUL);
AC_RM = quat2rotm(AC_EUL);
GH_RM = quat2rotm(GH_EUL);

SC_Q = rotm2eul(SC_RM,'YZX');
AC_Q = rotm2eul(AC_RM,'YZX');
GH_Q = rotm2eul(GH_RM,'YZY');

EL_x = angles(:,10);

res = [SC_Eul,AC_Eul,GH_Eul,EL_x];
end
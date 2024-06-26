function angles = rotxyz_sym(in1)
%ROTXYZ_SYM
%    ANGLES = ROTXYZ_SYM(IN1)

%    This function was generated by the Symbolic Math Toolbox version 9.3.
%    30-Jun-2024 17:51:02

q1 = in1(1,:);
q2 = in1(2,:);
q3 = in1(3,:);
q4 = in1(4,:);
q5 = in1(5,:);
q6 = in1(6,:);
q7 = in1(7,:);
q8 = in1(8,:);
q9 = in1(9,:);
t2 = cos(q1);
t3 = cos(q2);
t4 = cos(q3);
t5 = cos(q4);
t6 = cos(q5);
t7 = cos(q6);
t8 = cos(q7);
t9 = cos(q8);
t10 = cos(q9);
t11 = sin(q1);
t12 = sin(q2);
t13 = sin(q3);
t14 = sin(q4);
t15 = sin(q5);
t16 = sin(q6);
t17 = sin(q7);
t18 = sin(q8);
t19 = sin(q9);
t20 = t2.*t4;
t21 = t5.*t7;
t22 = t8.*t10;
t23 = t2.*t13;
t24 = t4.*t11;
t25 = t5.*t16;
t26 = t7.*t14;
t27 = t10.*t17;
t28 = t11.*t13;
t29 = t14.*t16;
t34 = t8.*t9.*t19;
t39 = t9.*t17.*t19;
t42 = t2.*t3.*t5.*t6;
t30 = t12.*t28;
t31 = t15.*t29;
t32 = t12.*t20;
t33 = t15.*t21;
t35 = t12.*t23;
t36 = t12.*t24;
t37 = t15.*t25;
t38 = t15.*t26;
t45 = -t39;
t46 = t27+t34;
t40 = -t30;
t41 = -t31;
t43 = -t32;
t44 = -t33;
t47 = t23+t36;
t48 = t24+t35;
t49 = t25+t38;
t50 = t26+t37;
t51 = t22+t45;
t52 = t20+t40;
t53 = t28+t43;
t54 = t21+t41;
t55 = t29+t44;
t56 = t2.*t3.*t50;
t57 = t6.*t14.*t48;
t64 = t48.*t49;
t58 = t15.*t53;
t59 = t6.*t7.*t53;
t60 = t2.*t3.*t55;
t61 = -t57;
t62 = t6.*t16.*t53;
t65 = t48.*t54;
t63 = -t62;
t66 = t42+t58+t61;
t68 = t59+t60+t64;
t67 = t46.*t66;
t69 = t56+t63+t65;
t70 = t18.*t19.*t68;
t71 = t51.*t69;
t72 = t67+t70+t71;
t73 = t72.^2;
t74 = -t73;
t75 = t74+1.0;
t76 = 1.0./sqrt(t75);
angles = [atan2(-t76.*(-t51.*(-t12.*t50+t3.*t13.*t54+t3.*t4.*t6.*t16)+t46.*(t3.*t4.*t15+t5.*t6.*t12+t3.*t6.*t13.*t14)+t18.*t19.*(t12.*t55-t3.*t13.*t49+t3.*t4.*t6.*t7)),-t76.*(t46.*(-t15.*t47+t6.*t14.*t52+t3.*t5.*t6.*t11)+t51.*(-t52.*t54+t3.*t11.*t50+t6.*t16.*t47)-t18.*t19.*(t49.*t52+t6.*t7.*t47-t3.*t11.*t55)));asin(t72);atan2(-t76.*(t9.*t68-t8.*t18.*t66+t17.*t18.*t69),t76.*(-t69.*(t8.*t19+t9.*t27)+t66.*(t9.*t22-t17.*t19)+t10.*t18.*t68))];
end

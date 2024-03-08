clc;clear all
data;
data_read;
t_end = 3.88;
initialConditions = zeros(1,18);

% Parametry pro kazde teleso
% c .. viskozita, k .. tuhost pruziny, I .. momenty setrvacnosti
% teles
% v simulinku je zadan pohyb pro jednotlive nezavisle souradnice q1..q9 a
% pomoci inverzni dynamiky dopocitany potrebne momentove sily
akt_GH = rand(1,13)*1;
akt_ST = rand(1,13)*1;

% OUT = sim('main',t_end);
OUT = sim('main_WU_motion.slx',t_end);
t = OUT.time.time;
fig = figure();
plot(t,OUT.q1(1,:),t,OUT.q2(1,:),t,OUT.q3(1,:),t,OUT.q10(1,:),'+',t,OUT.q11(1,:),'square',t,OUT.q12(1,:),'diamond')
title('UPPER BODY - MOMENTS')
xlabel('Time [s]')
legend('T1 SimScape [rad]','T2 SimScape [rad]','T3 SimScape [rad]','T1sim Sympy [rad]','T2sim Sympy [rad]','T3sim Sympy [rad]')
% 
fig = figure();
plot(t,OUT.q4(1,:),t,OUT.q5(1,:),t,OUT.q6(1,:),t,OUT.q13(1,:),'+',t,OUT.q14(1,:),'square',t,OUT.q15(1,:),'diamond')
title('MIDDLE BODY - MOMENTS')
xlabel('Time [s]')
legend('T4 SimScape [rad]','T5 SimScape [rad]','T6 SimScape [rad]','T4sim Sympy [rad]','T5sim Sympy [rad]','T6sim Sympy [rad]')

fig = figure();
plot(t,OUT.q7(1,:),t,OUT.q8(1,:),t,OUT.q9(1,:),t,OUT.q16(1,:),'+',t,OUT.q17(1,:),'square',t,OUT.q18(1,:),'diamond')
title('LOWER BODY - MOMENTS')
xlabel('Time [s]')
legend('T4 SimScape [rad]','T5 SimScape [rad]','T6 SimScape [rad]','T4sim Sympy [rad]','T5sim Sympy [rad]','T6sim Sympy [rad]')
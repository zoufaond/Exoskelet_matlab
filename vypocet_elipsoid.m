clear all

Tl = [-0.094965, 0.130202, -0.086922];
Te = [-0.05, 0.06, -0.03];

T_rel = Tl-Te;

h = 0.15

y_o = sqrt(T_rel(1)^2+T_rel(3)^2);
y_h = Tl(2)-Te(2);

r = sqrt((y_o^2))/(1-(y_h^2/h^2))
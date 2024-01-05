%% ellipsoid joint
% scapula frames initial position (measured)
S_00 = [-0.094965;+0.130202;-0.086922];
S_10 = [-0.062046;+0.163320;-0.091534];
S_20 = [-0.085749;+0.051924;-0.113461];
% initial position optimization
nvars = 2;
lb = [-pi;h_min];
ub = [+pi;h_max];

fun = @(X) dst(X(2),X(1),S_10,x_ej,y_ej,z_ej,psi_ej,mu_ej,phi_ej,r_ej,h_ej);
[X_10,fval_10] = ga(fun,nvars,[],[],[],[],lb,ub);
phi_10 = X_10(1);
h_10 = X_10(2);
r_10 = sqrt((r_ej^2)*(1-h_10^2/h_ej^2));

fun = @(X) dst(X(2),X(1),S_20,x_ej,y_ej,z_ej,psi_ej,mu_ej,phi_ej,r_ej,h_ej);
[X_20,fval_20] = ga(fun,nvars,[],[],[],[],lb,ub);
phi_20 = X_20(1);
h_20 = X_20(2);
r_20 = sqrt((r_ej^2)*(1-h_20^2/h_ej^2));

%% functions
% transformation matrices
function T = T_tx(x)
    T = [   1   0   0   x
            0   1   0   0
            0   0   1   0
            0   0   0   1   ];
end
function T = T_ty(y)
    T = [   1   0   0   0
            0   1   0   y
            0   0   1   0
            0   0   0   1   ];
end
function T = T_tz(z)
    T = [   1   0   0   0
            0   1   0   0
            0   0   1   z
            0   0   0   1   ];
end
function T = T_rx(phi_x)
    T = [   1           0           0           0
            0           +cos(phi_x) -sin(phi_x) 0
            0           +sin(phi_x) +cos(phi_x) 0
            0           0           0           1   ];
end
function T = T_ry(phi_y)
    T = [   +cos(phi_y) 0           +sin(phi_y) 0
            0           1           0           0
            -sin(phi_y) 0           +cos(phi_y) 0
            0           0           0           1   ];
end
function T = T_rz(phi_z)
    T = [   +cos(phi_z) -sin(phi_z) 0           0
            +sin(phi_z) +cos(phi_z) 0           0
            0           0           1           0
            0           0           0           1   ];
end
% ellipsoid joint
function d = dst(h,phi,S,x_ej,y_ej,z_ej,psi_ej,mu_ej,phi_ej,r_ej,h_ej)
    r = sqrt((r_ej^2)*(1-h^2/h_ej^2));
    T_t = T_tx(x_ej)*T_ty(y_ej)*T_tz(z_ej);
    T_r = T_ry(psi_ej)*T_rx(mu_ej)*T_ry(phi_ej);
    T_ej = T_rx(-pi/2)*T_rz(phi)*T_tz(h)*T_ry(pi/2)*T_tz(r);
    R = T_t*T_r*T_ej*[0;0;0;1];
    d = sqrt((S(1)-R(1))^2+(S(2)-R(2))^2+(S(3)-R(3))^2);
end

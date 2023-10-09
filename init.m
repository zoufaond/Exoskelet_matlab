%% ellipsoid joint
% scapula SHR frame initional position
% SHR_0 = [   -0.094965
%             +0.130202
%             -0.086922   ];

% initial position
h_0 = 0.5*h_ej;
r_0 = sqrt((r_ej^2)*(1-h_0^2/h_ej^2));
phi_0 = atan(-(r_ej*h_0)/(h_ej^2*sqrt(1-h_0^2/h_ej^2)));

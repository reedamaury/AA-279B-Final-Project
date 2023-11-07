function [J] = AA279jacobianshoot(odefun, y0in, tspan, dv)
% Returns 3x3 Jacobian for shooting method
%
% AA279 Function Library
% Last modified: 19 April 2018 by Andrew K. Barrows
%
% Jacobian describes sensitivity of final position to initial velocity
%
% State must be of the form [rx ry rz vx vy vz]'
%
% sample function call:
%     J = jacobianshoot(@earthorbitbody, [r1; v1], 75*60, 0.01)
%
% odefun  function handle for state derivative function
% y0in    state at which J should be evaluated
% tspan   trajectory start and end times [sec]
% dv      small velocity increment used for finite differences

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);

J = zeros(3, 3);

% Finite differences with +/- vx
y0 = y0in + [0; 0; 0; +dv; 0; 0];
[t,y] = ode113(odefun, tspan, y0, options);
rhigh = y(length(t),1:3)';
y0 = y0in + [0; 0; 0; -dv; 0; 0];
[t,y] = ode113(odefun, tspan, y0, options);
rlow = y(length(t),1:3)';
J(:, 1) = (rhigh-rlow)/(2*dv);

% Finite differences with +/- vy
y0 = y0in + [0; 0; 0; 0; +dv; 0];
[t,y] = ode113(odefun, tspan, y0, options);
rhigh = y(length(t),1:3)';
y0 = y0in + [0; 0; 0; 0; -dv; 0];
[t,y] = ode113(odefun, tspan, y0, options);
rlow = y(length(t),1:3)';
J(:, 2) = (rhigh-rlow)/(2*dv);

% Finite differences with +/- vz
y0 = y0in + [0; 0; 0; 0; 0; +dv];
[t,y] = ode113(odefun, tspan, y0, options);
rhigh = y(length(t),1:3)';
y0 = y0in + [0; 0; 0; 0; 0; -dv];
[t,y] = ode113(odefun, tspan, y0, options);
rlow = y(length(t),1:3)';
J(:, 3) = (rhigh-rlow)/(2*dv);

end % terminates MATLAB function





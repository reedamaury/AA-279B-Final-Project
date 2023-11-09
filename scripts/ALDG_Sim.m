clear; clc; close all;

% Define intial parameters for Starship HLS
starship_drymass = 85e5; % kg 
propellant_mass = 1e6; % kg 
m0 =  starship_drymass + propellant_mass; % kg
g = -1.62; % Moon gravity
Isp = 345;
vex = abs(Isp*g); 
Tmax = 2.2e6;
Tmin = Tmax*0.4;

% Define initial conditions (60 nautical miles above surface)
initial_conditions = [-20000; 0; 111120; 88.163; 0; -500; m0];
tspan = [0 500];  % Start and end times

% Define final conditions 
target_conditions = zeros(6,1);


options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9, 'Events', @eventsFcn);

% Solve the ODE system
[t, y, te, ye, ie] = ode113(@(t, y) Apollo_Lunar_Descent(t, y, tspan(2)-t, target_conditions,Tmax, Tmin, vex), tspan, initial_conditions,options);

% Extract the solution
rx = y(:, 1);
ry = y(:, 2);
rz = y(:, 3);
Vx = y(:, 4);
Vy = y(:, 5);
Vz = y(:, 6);
mass = y(:, 7);
rx(end-1)
rz(end-1)
figure;
plot3(rx, ry, rz);  % Plot with line and marker
xlabel('X position (m)');
ylabel('Y position (m)');
zlabel('Z position (m)');
title('Landing Vehicle Trajectory');
grid on;  % Add a grid
view(45, 45);  % Set a specific view angle

% Set axis limits
delta = 30;  % Adjust this value as needed
xlim([min(rx)-delta, max(rx)+delta]);
ylim([min(ry)-delta, max(ry)+delta]);
zlim([min(rz)-delta, max(rz)+delta]);
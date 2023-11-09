clear; clc;

% Define intial parameters 
starship_drymass = 85e5; % kg 
propellant_mass = 1e6; % kg 
m0 =  starship_drymass + propellant_mass; % kg
g = -1.62; % Moon gravity
Isp = 345;
vex = abs(Isp*g); 
Tmax = 2.2e6;
Tmin = Tmax*0.4;

% Define initial conditions 
initial_conditions = [-20000; 0; 111120; 88.171; 0; -500; m0];
tspan = [0 500];  % Start and end times

% Define final conditions 
target_conditions = zeros(6,1);

max_iterations = 2; % Arbitrary max number to ensure the loop doesn't run indefinitely
dv = 1e-4; % [m/s]
tolerance = 1; % 1 m
iteration = 0;

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9, 'Events', @eventsFcn);

% Solve the ODE system
[t, y, te, ye, ie] = ode113(@(t, y) Apollo_Lunar_Descent(t, y, tspan(2)-t, target_conditions,Tmax, Tmin, vex), tspan, initial_conditions,options);
rfinal_updated = y(end-1,1:3).'; % Final simulated position with updated v1
error_vector_updated = rfinal_updated - target_conditions(1:3);
error_magnitude = norm(error_vector_updated)
y0_updated = initial_conditions;

while (error_magnitude > tolerance || isnan(error_magnitude)) && iteration < max_iterations
       
    % Compute the Jacobian
    J = AA279jacobianshoot(y0_updated, tspan, dv, target_conditions, Tmax, Tmin, vex)
    
    % Compute the change in velocity using the Jacobian
    delta_v = J \ error_vector_updated;
    
    % Update v1
    v1_updated = y0_updated(4:6) + delta_v
    y0_updated(4:6) = v1_updated;

    % Run the simulation using the updated v1 value
    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9, 'Events', @eventsFcn);
    % Solve the ODE system
    [t, y_updated, te, ye, ie] = ode113(@(t, y) Apollo_Lunar_Descent(t, y, tspan(2)-t, target_conditions,Tmax, Tmin, vex), tspan, y0_updated, options);
    
    % Compute the new error vector
    rfinal_updated = y_updated(end-1,1:3)'; % Final position with updated v1
    error_vector_updated = rfinal_updated  - target_conditions(1:3);
    
    % Update the iteration count and error magnitude
    iteration = iteration + 1;
    error_magnitude = norm(error_vector_updated)
end

fprintf('Number of iterations required: %d\n', iteration);
fprintf('Final v1x: %f m/s\n', v1_updated(1));
fprintf('Final v1y: %f m/s\n', v1_updated(2));
fprintf('Final v1z: %f m/s\n', v1_updated(3));
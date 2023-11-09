clc ;
clear;

% terminal conditions/constraints
rfx= 0; % m
rfy=0; % m
rfz = 0; % m
vfx=0; % m/s
vfy=0; % m/s
vfz=0; % m/s

% initial conditions
rx(1) = -20000; % m 
ry(1) = 0; % m
rz(1) = 111120; % m
Vx(1) = 87.75; % m/s
Vy(1) = 0; % m/s
Vz(1) = -500; % m/s
starship_drymass = 85e5; % kg 
propellant_mass = 1e6; % kg 
m(1) =  starship_drymass + propellant_mass; % kg


% initial parameters
g = -1.62; % Moon gravity
Isp = 345;
vex = abs(Isp*g); 
Tmax = 2.2e6;
Tmin = Tmax*0.4;
tgo= 500;
h = 0.01;  % set the step size
x = [0:h:tgo];  % time interval
n = length(x)-1;

%%
%Function Definitions 
atx = @(~,atxcom)(atxcom);
aty = @(~,atycom)(atycom);
atz = @(~,atzcom)(atzcom);
m_dot = @(~, m, atcom)(-m.*atcom/vex);
rx_dot = @(~,Vx)(Vx); 
ry_dot = @(~, Vy)(Vy); 
rz_dot = @(~, Vz)(Vz); 


%%
%Apollo lunar descent guidance simulation
for i = 1:n
    atxcom(i) = (6*(rfx-rx(i)-((tgo)*Vx(i)))/(tgo).^2)-((2*(vfx-Vx(i)))/(tgo));
    atycom(i) = (6*(rfy-ry(i)-((tgo)*Vy(i)))/(tgo).^2)-((2*(vfy-Vy(i)))/(tgo));
    atzcom(i) = (6*(rfz-rz(i)-((tgo)*Vz(i)))/(tgo).^2)-((2*(vfz-Vz(i)))/(tgo))-g;
    atcom(i) = sqrt(atxcom(i).^2 + atycom(i).^2 + atzcom(i).^2);
    amin(i) = Tmin./m(i); 
    amax(i) = Tmax./m(i);
    
    %Physical Constraint Check (Maximum Allowable Acceleration/Thrust Given Engine Constraints)
    if  atcom(i) > amax(i)
        atxcom(i) = amax(i)*(atxcom(i)/atcom(i));
        atycom(i) = amax(i)*(atycom(i)/atcom(i));
        atzcom(i) = amax(i)*(atzcom(i)/atcom(i));
        atcom(i) = sqrt(atxcom(i).^2 + atycom(i).^2 + atzcom(i).^2);
        
        elseif atcom(i) < amin(i)
        atxcom(i) = amin(i)*(atxcom(i)/atcom(i));
        atycom(i) = amin(i)*(atycom(i)/atcom(i));
        atzcom(i) = amax(i)*(atzcom(i)/atcom(i));
        atcom(i) = sqrt(atxcom(i).^2 + atycom(i).^2 + atzcom(i).^2);
  
    end
    
    %4th Order Runge Kutta Numerical Integration
    k1 = atx(x(i),atxcom(i));
    k2 = atx(x(i)+.5*h,atxcom(i)+.5*k1*h);
    k3 = atx(x(i)+.5*h,atxcom(i)+.5*k2*h);
    k4 = atx(x(i)+h,atxcom(i)+k3*h);
    Vx(i+1) = Vx(i)+((k1+2*k2+2*k3+k4)/6)*h;  %+ normrnd(0,0.05);
    
    k1 = aty(x(i), atycom(i));
    k2 = aty(x(i)+.5*h,atycom(i)+.5*k1*h);
    k3 = aty(x(i)+.5*h,atycom(i)+.5*k2*h);
    k4 = aty(x(i)+h,atycom(i)+k3*h);
    Vy(i+1) = Vy(i)+((k1+2*k2+2*k3+k4)/6)*h ; % + normrnd(0,0.05);

    k1 = atz(x(i),atzcom(i));
    k2 = atz(x(i)+.5*h,atzcom(i)+.5*k1*h);
    k3 = atz(x(i)+.5*h,atzcom(i)+.5*k2*h);
    k4 = atz(x(i)+h,atzcom(i)+k3*h);
    Vz(i+1) = Vz(i)+((k1+2*k2+2*k3+k4)/6)*h;  %+ normrnd(0,0.05);
      
    
    k1 =  m_dot(x(i),m(i),atcom(i));
    k2 =  m_dot(x(i)+.5*h,m(i)+.5*k1*h, atcom(i));
    k3 =  m_dot(x(i)+.5*h,m(i)+.5*k2*h, atcom(i));
    k4 =  m_dot(x(i)+h,m(i)+k3*h, atcom(i));
    m(i+1) = m(i)+((k1+2*k2+2*k3+k4)/6)*h;

    k1 = rx_dot(x(i),Vx(i));
    k2 = rx_dot(x(i)+.5*h,Vx(i)+.5*k1*h);
    k3 = rx_dot(x(i)+.5*h,Vx(i)+.5*k2*h);
    k4 = rx_dot(x(i)+h,Vx(i)+k3*h);
    rx(i+1) = rx(i)+((k1+2*k2+2*k3+k4)/6)*h; % + normrnd(0,0.5);


    k1 = ry_dot(x(i),Vy(i));
    k2 = ry_dot(x(i)+.5*h,Vy(i)+.5*k1*h);
    k3 = ry_dot(x(i)+.5*h,Vy(i)+.5*k2*h);
    k4 = ry_dot(x(i)+h,Vy(i)+k3*h);
    ry(i+1) = ry(i)+((k1+2*k2+2*k3+k4)/6)*h; % + normrnd(0,0.5);
    

    k1 = rz_dot(x(i),Vz(i));
    k2 = rz_dot(x(i)+.5*h,Vz(i)+.5*k1*h);
    k3 = rz_dot(x(i)+.5*h,Vz(i)+.5*k2*h);
    k4 = rz_dot(x(i)+h,Vz(i)+k3*h);
    rz(i+1) = rz(i)+((k1+2*k2+2*k3+k4)/6)*h; % + normrnd(0,0.5);
    
    %Ensure Altitude is Greater Than 0
    if rz(i+1) < 0 
        break 
    end

    tgo = tgo-h;
    if tgo < 0 
        break
    end
    
end


Pc = m(1)-m(end-1)
propellant_mass
v = sqrt(Vx.^2+Vy.^2 + Vz.^2);
T = m(1:end-1).*atcom(1:end);
rx(end)
% %% Plots
% figure(1)
% plot(x(1:length(ry)-1),ry(1:end-1))
% title('Altitude vs. Time')
% xlabel('Time (s)') 
% ylabel('y(t) (m)')
% xlim([0 110])
% ylim([0 4500])


% figure(2)
% plot(x(1:length(v)-1),v(1:end-1))
% title('Velocity Magnitude vs. Time')
% xlabel('Time (s)')
% ylabel('Velocity (m/s)')
% xlim([0 110]) 
% ylim([0 600])

% figure(2)
% plot(x(1:length(T)),m(1:end-1))
% hold on
% title('Mass Profile')
% xlabel('Time (s)')
% ylabel('Mass (kg)')
% xlim([0 50])
% ylim([1450 2000])
% hold off



figure;
plot3(rx, ry, rz, '-o');  % Plot with line and marker
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



% figure; 
% scatter3(rx, rz, ry)
% figure(4)
% plot(x(1:length(Vx)-1), Vx(1:end-1))
% hold on
% plot(x(1:length(Vy)-1), Vy(1:end-1))
% title('X and Y Velocity vs. Time')
% xlabel('Time (s)')
% ylabel('Velocity (m/s)')
% legend('X Velocity (m/s)', 'Y Velocity (m/s)')
% hold off
% 
% figure(5)
% plot(x(1:length(T)),T)
% hold on
% title('Thrust Profile')
% xlabel('Time (s)')
% ylabel('Thrust (N)')
% xlim([0 50])
% ylim([0 18600])
% hold off

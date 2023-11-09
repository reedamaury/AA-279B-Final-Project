function dydt = Apollo_Lunar_Descent(t, y, tgo, target_conditions, Tmax, Tmin, vex)
    mu_moon = 4903;  % Moon's gravitational parameter [km^3/s^2]
    r_rel = y(1:3)./1000; 
    r_rel_mag = norm(r_rel); 
    rfx = target_conditions(1);
    rfy = target_conditions(2);
    rfz = target_conditions(3);
    vfx = target_conditions(4);
    vfy = target_conditions(5);
    vfz = target_conditions(6);

    atxcom = (6*(rfx-y(1)-((tgo)*y(4)))/(tgo).^2)-((2*(vfx-y(4)))/(tgo)); 
    atycom = (6*(rfy-y(2)-((tgo)*y(5)))/(tgo).^2)-((2*(vfy-y(5)))/(tgo));
    atzcom = (6*(rfz-y(3)-((tgo)*y(6)))/(tgo).^2)-((2*(vfz-y(6)))/(tgo)) + mu_moon/r_rel_mag^3 * r_rel(3);
    atcom = sqrt(atxcom.^2 + atycom.^2 + atzcom.^2);
    amin = Tmin./y(7); 
    amax = Tmax./y(7);
    
    % Physical Constraint Check (Maximum Allowable Acceleration/Thrust Given Engine Constraints)
    if  atcom > amax
        atxcom = amax*(atxcom/atcom);
        atycom = amax*(atycom/atcom);
        atzcom = amax*(atzcom/atcom);
        atcom = sqrt(atxcom.^2 + atycom.^2 + atzcom.^2);
        
        elseif atcom < amin
        atxcom = amin*(atxcom/atcom);
        atycom = amin*(atycom/atcom);
        atzcom = amax*(atzcom/atcom);
        atcom = sqrt(atxcom.^2 + atycom.^2 + atzcom.^2);
  
    end
    
    dydt = zeros(7,1);  % initialize the output
    dydt(1:3) = y(4:6); % velocities
    dydt(4) = atxcom;
    dydt(5) = atycom;
    dydt(6) = atzcom; 
    dydt(7) = (-y(7).*atcom/vex);

end
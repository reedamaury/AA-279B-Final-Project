function [r_eci, v_eci] = launcherPositionVelocity( dayOfYear, hour, min, sec)
    % Given constants
    muEarth = 398600.4418; % [km3/sec2]
    rEarth = 6378.137; % [km]
    omegaEarth = 0.0000729211585530; % [rad/sec]
    theta_g0 = 6.7045617 * (2 * pi) / 24; % [rad] Initial Greenwich sidereal time converted from hours to radians
    e = 0.081819221456; % eccentricity of Earth (oblate ellipsoid)

    % Given launcher data (Kennedy Space Center)
    longitude = -80.651070 * (pi / 180); % [rad]
    latitude = 28.573469 * (pi / 180); % [rad]
    elevation = 0.003; % [km]

    % Calculate time from the beginning of 2014 in seconds
    sec_in_day = 86164;
    elapsed_time = (dayOfYear - 1) * sec_in_day + hour * 3600 + min * 60 + sec;

    % Compute GST for the desired time
    theta_g = theta_g0 + omegaEarth * elapsed_time;
    theta_g = mod(theta_g, 2 * pi); % Ensure 0 <= theta_g < 2*pi

    % Compute launcher position in ECEF (rotating earth)
    N = rEarth / sqrt(1 - e^2 * sin(latitude)^2); % Earth ellipsoid model, WGS 84
    x_ecef = (N + elevation) * cos(latitude) * cos(longitude);
    y_ecef = (N + elevation) * cos(latitude) * sin(longitude);
    z_ecef = (N * (1 - e^2) + elevation) * sin(latitude);
    r_ecef = [x_ecef; y_ecef; z_ecef];

    % Rotate from ECEF to ECI
    R = [cos(theta_g), sin(theta_g), 0;
        -sin(theta_g), cos(theta_g), 0;
        0, 0, 1]';
    r_eci = R * r_ecef;

    % Compute inertial velocity of the launcher in ECI
    v_ecef = [0; 0; 0]; 
    v_eci = R' * (v_ecef + cross([0; 0; -omegaEarth], r_ecef));

end


function [value, isterminal, direction] = eventsFcn(t, y)
    value = y(3);      % We're interested in the height above the surface, y(3)
    isterminal = 1;    % Stop the integration when value crosses zero
    direction = -1;    % The zero is approached from above (height decreasing)
end

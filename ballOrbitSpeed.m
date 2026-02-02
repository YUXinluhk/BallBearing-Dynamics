function omega_m = ballOrbitSpeed(omega, sin_alp_i, cos_alp_i,...
    sin_alp_o, cos_alp_o, gamma_tick, rc_id)
%Calculates the orbital speed of the balls.
switch rc_id
    case 0 %inner race control
        omega_m=omega*(((cos_alp_i.*cos_alp_o+sin_alp_i.*sin_alp_o)...
            -gamma_tick.*cos_alp_o)./...
            (1+(cos_alp_i.*cos_alp_o+sin_alp_i.*sin_alp_o)));
    case 1 %outer race control
        omega_m=omega*((1-gamma_tick.*cos_alp_i)./...
            (1+(cos_alp_i.*cos_alp_o+sin_alp_i.*sin_alp_o)));
end
%trigonometric identity: cosd(alpha_i-alpha_o) =
%                  cosd(alpha_i)*cosd(alpha_o)+sind(alpha_i)*sind(alpha_o)
end

function M_g = gyroMoment(obj, J, omega, sin_alp_i, cos_alp_i,...
    sin_alp_o, cos_alp_o, gamma_tick, rc_id)
switch rc_id
    case 0 %inner race control
        beta=atan2d(sin_alp_i,cos_alp_i-gamma_tick);
    case 1 %outer race control
        beta=atan2d(sin_alp_o,cos_alp_o+gamma_tick);
end
omega_m =  ballOrbitSpeed(omega, sin_alp_i, cos_alp_i,...
    sin_alp_o, cos_alp_o, gamma_tick, rc_id);
omega_R = ballSpinSpeed(sin_alp_i, cos_alp_i, sin_alp_o,...
    cos_alp_o, omega, gamma_tick, beta);
M_g=J*omega_R.*omega_m.*sind(beta);
end

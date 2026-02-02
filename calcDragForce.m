function F_drag = calcDragForce(omega_m, D, dm, lub)
%calcDragForce Calculates aerodynamic drag force on the ball
%   Based on Equation 2-32: Fd = 1/32 * pi * Cd * rho * (D * dm * omega_m)^2
%   Inputs:
%       omega_m: Orbital speed [rad/s]
%       D: Ball Diameter [m]
%       dm: Pitch Diameter [m]
%       lub: Lubricant struct (rho)
%   Outputs:
%       F_drag: Drag Force [N] (acting opposite to orbital motion)

Cd = 1.0; % Drag coefficient (approximate for sphere)
% Note: Paper Eq 2-32 uses Cd and rho_eff (mixture density)
% We rely on lub.rho which implies effective density if pre-calculated,
% or we should use X (grease fraction). For now assuming lub.rho is effective.

rho_eff = lub.rho;

% Eq 2-32
% Fd = (pi * Cd * rho * (D * dm * omega_m)^2) / 32
% Note: Velocity of ball center V = dm/2 * omega_m
% Area = pi * (D/2)^2 ??
% Standard Drag: 0.5 * rho * V^2 * Cd * A
% V = 0.5 * dm * wm
% A = pi * D^2 / 4
% F = 0.5 * rho * (0.25*dm^2*wm^2) * Cd * (0.25*pi*D^2)
%   = (1/32) * rho * Cd * pi * D^2 * dm^2 * wm^2
% Matches the formula structure exactly.

F_drag = (pi * Cd * rho_eff * (D * dm * omega_m)^2) / 32;

end

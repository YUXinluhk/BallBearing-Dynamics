function [aftStruct,rc_inner] = addSigmaMaxAlpha(obj, befStruct, geo, phys,...
    load)
aftStruct=befStruct;
A_1 = raceCenterDistAx(aftStruct.delta_a,...
    geo, load);
A_2 = raceCenterDistRad(aftStruct.delta_ry,...
    aftStruct.delta_rz, geo);

[cos_alp_i, cos_alp_o, sin_alp_i, sin_alp_o] =...
    trigSincosd(obj, A_1, A_2,...
    aftStruct.X_1, aftStruct.X_2,...
    aftStruct.delta_i, aftStruct.delta_o,...
    geo.D, geo.f_o, geo.f_i);

rho_ci = curvInner(cos_alp_i, geo.f_i, geo.D,...
    geo.d_m);
rho_co = curvOuter(cos_alp_o, geo.f_o, geo.D,...
    geo.d_m);

[Q_i, aftStruct.sigma_max_i, a_i, E_i] = calcStress(obj, ...
    aftStruct.delta_i, rho_ci, phys);
[Q_o, aftStruct.sigma_max_o, a_o, E_o] = calcStress(obj, ...
    aftStruct.delta_o, rho_co, phys);

aftStruct.alpha_o=atan2d(sin_alp_o,cos_alp_o);
aftStruct.alpha_i=atan2d(sin_alp_i,cos_alp_i);

aftStruct.Q_i=Q_i;
aftStruct.Q_o=Q_o;
rc_trig=cos_alp_i.*cos_alp_o+sin_alp_i.*sin_alp_i;

rc_inner=rcCheck(Q_o, a_o, E_o, rc_trig, Q_i, a_i, E_i);

% Quasi-Dynamic: Beta and Omega are solved variables. Do not overwrite with Kinematic approximation.
% omega=(load.n/60)*2*pi;
% rc_id=not(rc_inner);
% switch rc_id
%     case 0 %inner race control
%         beta=atan2d(sin_alp_i,cos_alp_i-geo.gamma_tick);
%     case 1 %outer race control
%         beta=atan2d(sin_alp_o,cos_alp_o+geo.gamma_tick);
% end
% aftStruct.omega_m = ballOrbitSpeed(omega, sin_alp_i,...
%     cos_alp_i, sin_alp_o, cos_alp_o, geo.gamma_tick, rc_id);
% aftStruct.omega_R = ballSpinSpeed(sin_alp_i, cos_alp_i,...
%     sin_alp_o, cos_alp_o, omega, geo.gamma_tick, beta);
% aftStruct.beta = beta;

% Force rc_inner to false to prevent re-running the solver with different control assumption
rc_inner = 0;

%This is good for debugging:
%clipboard('copy', [cos_alp_o(1) cos_alp_i(1);sin_alp_o(1) sin_alp_i(1); a_o(1) a_i(1)]);
end

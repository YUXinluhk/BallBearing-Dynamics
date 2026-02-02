function F_c = centrifugalForce(obj, d_m, omega, sin_alp_i,...
    cos_alp_i, sin_alp_o, cos_alp_o, gamma_tick, m_ball, rc_id)
%centrifugalForce 计算球的离心力
%   根据论文公式 2-49: F_η = 0.5 * m * ω_m² * d_m
%   其中 ω_m 是球的公转角速度, d_m 是节圆直径
%   物理意义：离心力 = m * r * ω² = m * (d_m/2) * ω_m²
%                    = 0.5 * m * d_m * ω_m²
omega_m = ballOrbitSpeed(omega, sin_alp_i, cos_alp_i,...
    sin_alp_o, cos_alp_o, gamma_tick, rc_id);
% 论文公式 2-49: F = 0.5 * m * ω_m² * d_m
F_c = 0.5 * m_ball * d_m * omega_m.^2;
end

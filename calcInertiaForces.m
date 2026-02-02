function [F_inertia, M_gyro] = calcInertiaForces(j, phys, geo, omega_m_all, omega_vec_all, delta_psi)
%calcInertiaForces 计算球的惯性力和陀螺力矩
%   完整实现论文公式 2-49, 2-50, 2-51
%
%   输入:
%       j: 球索引 (1 to z)
%       phys: 物理参数结构体 (含 m_ball, J_ball)
%       geo: 几何参数结构体 (含 d_m, z)
%       omega_m_all: 所有球的公转速度向量 [rad/s]
%       omega_vec_all: 所有球的自转分量矩阵 [z x 3] (omega_x, omega_y, omega_z)
%       delta_psi: 球间角间距 [rad] (默认 2π/z)
%
%   输出:
%       F_inertia: 惯性力 [F_b, F_eta, F_tau] [N] (公式 2-49)
%       M_gyro: 陀螺力矩 [M_x, M_y, M_z] [N.m] (公式 2-50)

if nargin < 6
    delta_psi = 2 * pi / geo.z;  % 默认均匀分布
end

% 提取参数
m = phys.m_ball;                % 球质量 [kg]
I = phys.J_ball;                % 球转动惯量 [kg.m^2] (I = mD²/10)
d_m = geo.d_m;                  % 节圆直径 [m]
z = geo.z;                      % 球数

% 当前球的公转速度
omega_m_j = omega_m_all(j);

% 当前球的自转分量
omega_x_j = omega_vec_all(j, 1);
omega_y_j = omega_vec_all(j, 2);
omega_z_j = omega_vec_all(j, 3);

% 获取相邻球索引 (循环边界)
j_prev = mod(j - 2, z) + 1;  % j-1 (循环)
j_next = mod(j, z) + 1;      % j+1 (循环)

% ========== 公式 2-49: 惯性力 ==========
%
% F_b = 0 (副法线方向，恒速时为零)
% F_η = 0.5 * m * ω_mj² * d_m (离心力，主法线方向)
% F_τ = 0.5 * m * ω_mj * d_m * (ω_m(j+1) - ω_m(j-1)) / ΔΨ (切向加速力)

F_b = 0;

F_eta = 0.5 * m * omega_m_j^2 * d_m;

% 切向加速力 (非稳态项)
omega_m_next = omega_m_all(j_next);
omega_m_prev = omega_m_all(j_prev);
if abs(delta_psi) > 1e-12
    F_tau = 0.5 * m * omega_m_j * d_m * (omega_m_next - omega_m_prev) / delta_psi;
else
    F_tau = 0;
end

F_inertia = [F_b; F_eta; F_tau];

% ========== 公式 2-51: 角速度微分 ==========
%
% ω̇_xj ≈ 0.5 * ω_mj * (ω_x(j+1) - ω_x(j-1)) / ΔΨ
% ω̇_yj ≈ 0.5 * ω_mj * (ω_y(j+1) - ω_y(j-1)) / ΔΨ
% ω̇_zj ≈ 0.5 * ω_mj * (ω_z(j+1) - ω_z(j-1)) / ΔΨ

omega_x_next = omega_vec_all(j_next, 1);
omega_x_prev = omega_vec_all(j_prev, 1);
omega_y_next = omega_vec_all(j_next, 2);
omega_y_prev = omega_vec_all(j_prev, 2);
omega_z_next = omega_vec_all(j_next, 3);
omega_z_prev = omega_vec_all(j_prev, 3);

if abs(delta_psi) > 1e-12
    omega_dot_x = 0.5 * omega_m_j * (omega_x_next - omega_x_prev) / delta_psi;
    omega_dot_y = 0.5 * omega_m_j * (omega_y_next - omega_y_prev) / delta_psi;
    omega_dot_z = 0.5 * omega_m_j * (omega_z_next - omega_z_prev) / delta_psi;
else
    omega_dot_x = 0;
    omega_dot_y = 0;
    omega_dot_z = 0;
end

% ========== 公式 2-50: 陀螺力矩 (欧拉方程) ==========
%
% M_xj = -I * ω̇_xj
% M_yj = -I * (ω̇_yj + ω_mj * ω_zj)
% M_zj = -I * (ω̇_zj - ω_mj * ω_yj)

M_x = -I * omega_dot_x;
M_y = -I * (omega_dot_y + omega_m_j * omega_z_j);
M_z = -I * (omega_dot_z - omega_m_j * omega_y_j);

M_gyro = [M_x; M_y; M_z];

end

function [F_cy, F_cz, M_cx] = calcCageHydrodynamics(cage, lub, geo)
%calcCageHydrodynamics 计算保持架-引导套圈流体动力学
%   实现论文公式 2-25 ~ 2-31
%
%   输入:
%       cage: 保持架结构体 (含 omega_c, epsilon_c, phi_c 等)
%       lub: 润滑参数结构体
%       geo: 几何参数结构体
%
%   输出:
%       F_cy, F_cz: 保持架所受流体动压力 [N] (全局坐标)
%       M_cx: 剪切摩擦力矩 [N.m]

% 提取保持架参数
omega_c = cage.omega_c;             % 保持架转速 [rad/s]
R1 = cage.R_outer;                  % 保持架外圆表面半径 [m]
L = cage.guide_width;               % 保持架引导面宽度 [m]
C1 = cage.guide_clearance;          % 保持架引导间隙 [m]
epsilon_c = cage.epsilon_c;         % 相对偏心率 [-] (公式 2-27)
delta_yc = cage.delta_yc;           % 保持架Y方向偏移 [m]
delta_zc = cage.delta_zc;           % 保持架Z方向偏移 [m]

% 提取润滑参数
eta0 = lub.eta_0;                   % 基础粘度 [Pa.s]

% ========== 公式 2-26: 保持架拖动速度 ==========
u1 = R1 * omega_c;

% ========== 公式 2-29: 保持架滑移速度 ==========
V1 = R1 * omega_c;

% ========== 公式 2-31: 最大油膜厚度角度 ==========
if abs(delta_yc) < 1e-12
    phi_c = pi/2 * sign(delta_zc);
else
    phi_c = atan2(delta_zc, delta_yc);
end

% ========== 公式 2-25: 流体动压力 (保持架局部坐标) ==========
% F'_cy = ± η₀·u₁·L³·ε² / [C₁²·(1-ε²)²]
% F'_cz = ∓ π·η₀·u₁·L³·ε / [4·C₁²·(1-ε²)^1.5]

% 防止除零
epsilon_safe = min(abs(epsilon_c), 0.99);
denom1 = C1^2 * (1 - epsilon_safe^2)^2;
denom2 = 4 * C1^2 * (1 - epsilon_safe^2)^1.5;

if denom1 < 1e-18 || denom2 < 1e-18
    F_prime_cy = 0;
    F_prime_cz = 0;
else
    F_prime_cy = eta0 * u1 * L^3 * epsilon_safe^2 / denom1;
    F_prime_cz = -pi * eta0 * u1 * L^3 * epsilon_safe / denom2;
end

% ========== 公式 2-28: 剪切摩擦力矩 ==========
% M'_cx = 2π·η₀·V₁·R₁·L / [C₁·√(1-ε²)]

denom3 = C1 * sqrt(1 - epsilon_safe^2);
if abs(denom3) < 1e-18
    M_prime_cx = 0;
else
    M_prime_cx = 2 * pi * eta0 * V1 * R1 * L / denom3;
end

% ========== 公式 2-30: 坐标变换到全局坐标系 ==========
% [M_cx]   [1      0        0     ] [M'_cx]
% [F_cy] = [0   cos(φ_c) -sin(φ_c)] [F'_cy]
% [F_cz]   [0   sin(φ_c)  cos(φ_c)] [F'_cz]

M_cx = M_prime_cx;
F_cy = cos(phi_c) * F_prime_cy - sin(phi_c) * F_prime_cz;
F_cz = sin(phi_c) * F_prime_cy + cos(phi_c) * F_prime_cz;

end

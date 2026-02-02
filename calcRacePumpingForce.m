function [F_R, F_S, F_H] = calcRacePumpingForce(lub, geo, u_slip, h_c, R_xi, R_eta, r_s, contact_type)
%calcRacePumpingForce 计算球-滚道入口区泵吸作用力
%   实现论文公式 2-41 ~ 2-48
%
%   输入:
%       lub: 润滑参数结构体
%       geo: 几何参数结构体
%       u_slip: 滑动速度 [u_xi, u_eta] [m/s]
%       h_c: 油膜厚度 [m]
%       R_xi, R_eta: 等效曲率半径 [m]
%       r_s: 入口区半径 [m] (通常取接触椭圆半轴)
%       contact_type: 'inner' 或 'outer'
%
%   输出:
%       F_R: 滚动泵吸力 [F_R_xi, F_R_eta] [N] (公式 2-41)
%       F_S: 滑动泵吸力 [F_S_xi, F_S_eta] [N] (公式 2-42)
%       F_H: 流体动压力 [F_H_xi, F_H_eta] [N] (公式 2-43)

if nargin < 8
    contact_type = 'inner';
end

% 提取参数
eta0 = lub.eta_0;
D = geo.D;

% 防止除零
R_xi = max(R_xi, 1e-9);
R_eta = max(R_eta, 1e-9);
h_c = max(h_c, 1e-9);
r_s = max(r_s, 1e-9);

% 椭圆比率
k = R_xi / R_eta;
k = max(k, 1e-3);

% 滑动速度分量
u_xi = u_slip(1);
u_eta = u_slip(2);

% 防止除零
if abs(u_xi) < 1e-12
    u_xi = 1e-12 * sign(u_xi + 1e-15);
end

% ========== 公式 2-44: C_0 系数 ==========
% C_0(i/o)j = η₀·u_ξj·√(R_ξ·R_η) × √[(3+2k)⁻² + (u_η/u_ξ)²·(3+2k⁻¹)⁻²·k⁻¹]

term1 = (3 + 2*k)^(-2);
term2 = (u_eta / u_xi)^2 * (3 + 2/k)^(-2) / k;

C_0 = eta0 * abs(u_xi) * sqrt(R_xi * R_eta) * sqrt(term1 + term2);

% ========== 公式 2-45: θ 角度 ==========
% θ(i/o)j = arctan[(3+2k) / (k^0.5·(3+2k⁻¹)) · (u_η/u_ξ)]

theta = atan((3 + 2*k) / (sqrt(k) * (3 + 2/k)) * (u_eta / u_xi));

% ========== 公式 2-46: t 参数 ==========
% t(i/o) = r_s · √[(cos²θ + sin²θ/k) / (2h·R_ξ)]

cos_theta = cos(theta);
sin_theta = sin(theta);

t_inner = (cos_theta^2 + sin_theta^2 / k) / (2 * h_c * R_xi);
t = r_s * sqrt(max(t_inner, 0));

% 防止 log(0)
t = max(t, 1e-6);

% ========== 公式 2-47: F̄_R 滚动力系数 ==========
% F̄_R(i/o)j = { 28.59·ln(t) - 10.10,  t ≤ 5
%             { 36.57·ln(t) - 22.85,  t > 5

if t <= 5
    F_bar_R = 28.59 * log(t) - 10.10;
else
    F_bar_R = 36.57 * log(t) - 22.85;
end
F_bar_R = max(F_bar_R, 0);  % 确保非负

% ========== 公式 2-48: F̄_S 滑动力系数 ==========
% F̄_S(i/o)j ≈ 0

F_bar_S = 0;

% ========== 公式 2-41: 滚动泵吸力 ==========
% F_R_ξ(i/o) = 0.5·C_0(i/o)·F̄_R(i/o)·cos(θ_i/o)
% F_R_η(i/o) = 0.5·C_0(i/o)·F̄_R(i/o)·√(R_ξ·R_η)·sin(θ_i/o)

F_R_xi = 0.5 * C_0 * F_bar_R * cos_theta;
F_R_eta = 0.5 * C_0 * F_bar_R * sqrt(R_xi * R_eta) * sin_theta;

F_R = [F_R_xi; F_R_eta];

% ========== 公式 2-42: 滑动泵吸力 ==========
% F_S_ξ(i/o) = F̄_S(i/o)·η₀·u_s(i/o)ξ·√(R_ξ·R_η)
% F_S_η(i/o) = F̄_S(i/o)·η₀·u_s(i/o)η·√(R_ξ·R_η)

sqrt_R = sqrt(R_xi * R_eta);
F_S_xi = F_bar_S * eta0 * u_xi * sqrt_R;
F_S_eta = F_bar_S * eta0 * u_eta * sqrt_R;

F_S = [F_S_xi; F_S_eta];

% ========== 公式 2-43: 流体动压力 (作用于球心) ==========
% F_H(i/o)_ξ = 2·C_0(i/o)·F̄_R(i/o)·R_ξ(i/o)·cos(θ_i/o) / D
% F_H(i/o)_η = 2·C_0(i/o)·F̄_R(i/o)·R_η(i/o)·√(R_ξ·R_η)·sin(θ_i/o) / D

F_H_xi = 2 * C_0 * F_bar_R * R_xi * cos_theta / D;
F_H_eta = 2 * C_0 * F_bar_R * R_eta * sqrt_R * sin_theta / D;

F_H = [F_H_xi; F_H_eta];

end

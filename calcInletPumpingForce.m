function [P_R, P_S, rho_eff] = calcInletPumpingForce(lub, geo, omega_m, u_slip, h_c, contact_type)
%calcInletPumpingForce 计算接触面入口区泵吸作用力
%   实现论文公式 2-33 ~ 2-40
%
%   输入:
%       lub: 润滑参数结构体
%       geo: 几何参数结构体
%       omega_m: 公转角速度 [rad/s]
%       u_slip: 滑动速度 [u_xi, u_eta] [m/s]
%       h_c: 油膜厚度 [m]
%       contact_type: 'ball-cage' 或 'ball-race'
%
%   输出:
%       P_R: 滚动泵吸力 [P_R_xi, P_R_eta] [N]
%       P_S: 滑动泵吸力 [P_S_xi, P_S_eta] [N]
%       rho_eff: 等效密度 [kg/m^3]

if nargin < 6
    contact_type = 'ball-race';
end

% ========== 公式 2-33: 等效密度 ==========
% ρ_e = X * ρ_grease + (1-X) * ρ_air
%
% 其中 X 为润滑脂比例系数

if isfield(lub, 'grease_ratio')
    X = lub.grease_ratio;
else
    X = 0.3;  % 默认 30% 润滑脂
end

rho_grease = lub.rho;       % 润滑脂密度 [kg/m^3]
rho_air = 1.2;              % 空气密度 [kg/m^3]

rho_eff = X * rho_grease + (1 - X) * rho_air;

% 提取润滑参数
eta0 = lub.eta_0;           % 基础粘度 [Pa.s]

% 根据接触类型选择参数
switch contact_type
    case 'ball-cage'
        % 球-保持架接触
        R_xi = geo.D / 2;                   % 短轴曲率半径
        R_eta = geo.D_p / 2;                % 长轴曲率半径 (兜孔)
        D_c = geo.d_m + geo.D;              % 外径
        d_c = geo.d_m - geo.D;              % 内径
    case 'ball-race'
        % 球-滚道接触
        R_xi = geo.D / 2;                   % 短轴曲率半径
        R_eta = geo.f_i * geo.D;            % 长轴曲率半径 (滚道曲率)
        D_c = geo.d_m;
        d_c = geo.d_m;
    otherwise
        R_xi = geo.D / 2;
        R_eta = geo.D / 2;
        D_c = geo.d_m;
        d_c = geo.d_m;
end

% 防止除零
R_xi = max(R_xi, 1e-9);
R_eta = max(R_eta, 1e-9);
h_c = max(h_c, 1e-9);

% 椭圆比率
k_p = R_xi / R_eta;
k_p = max(k_p, 1e-3);

% 滑动速度分量
u_xi = u_slip(1);
u_eta = u_slip(2);

% 防止除零
if abs(u_xi) < 1e-12
    u_xi = 1e-12 * sign(u_xi + 1e-15);
end

% ========== 公式 2-36: C_0p 系数 ==========
% C_0pj = η₀·u_ξj·√(R_pξ·R_pη) × √[(3+2k)⁻² + (u_η/u_ξ)²·(3+2k⁻¹)⁻²·k⁻¹]

term1 = (3 + 2*k_p)^(-2);
term2 = (u_eta / u_xi)^2 * (3 + 2/k_p)^(-2) / k_p;

C_0p = eta0 * abs(u_xi) * sqrt(R_xi * R_eta) * sqrt(term1 + term2);

% ========== 公式 2-37: θ_p 角度 ==========
% θ_pj = arctan[(3+2k) / (k^0.5·(3+2k⁻¹)) · (u_η/u_ξ)]

theta_p = atan((3 + 2*k_p) / (sqrt(k_p) * (3 + 2/k_p)) * (u_eta / u_xi));

% ========== 公式 2-38: t_p 参数 ==========
% t_p = 0.25(D_c - d_c) × √[(cos²θ + sin²θ/k) / (2h_c·R_ξ)]

cos_theta = cos(theta_p);
sin_theta = sin(theta_p);

t_p_inner = (cos_theta^2 + sin_theta^2 / k_p) / (2 * h_c * R_xi);
t_p = 0.25 * abs(D_c - d_c) * sqrt(max(t_p_inner, 0));

% 防止 log(0)
t_p = max(t_p, 1e-6);

% ========== 公式 2-39: P̄_R 滚动压力参数 ==========
% P̄_Rj = 34.74·ln(t_p) - 27.6

P_bar_R = 34.74 * log(t_p) - 27.6;
P_bar_R = max(P_bar_R, 0);  % 确保非负

% ========== 公式 2-40: P̄_S 滑动压力参数 ==========
% P̄_Sj = 0.26·P̄_Rj + 10.90

P_bar_S = 0.26 * P_bar_R + 10.90;

% ========== 公式 2-34: 滚动泵吸力 ==========
% P_R_ξ = 0.5·C_0p·P̄_R·cos(θ_p)
% P_R_η = 0.5·C_0p·P̄_R·√(R_ξ/R_η)·sin(θ_p)

P_R_xi = 0.5 * C_0p * P_bar_R * cos_theta;
P_R_eta = 0.5 * C_0p * P_bar_R * sqrt(R_xi / R_eta) * sin_theta;

P_R = [P_R_xi; P_R_eta];

% ========== 公式 2-35: 滑动泵吸力 ==========
% P_S_ξ = P̄_S·η₀·u_ξ·√(R_ξ·R_η)
% P_S_η = P̄_S·η₀·u_η·√(R_ξ·R_η)

sqrt_R = sqrt(R_xi * R_eta);
P_S_xi = P_bar_S * eta0 * u_xi * sqrt_R;
P_S_eta = P_bar_S * eta0 * u_eta * sqrt_R;

P_S = [P_S_xi; P_S_eta];

end

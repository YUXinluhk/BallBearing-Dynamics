function [v_slide_i, v_slide_o, omega_s_i, omega_s_o, omega_vec, kin_params] = calcBallKinematics(...
    omega_m, omega_b, beta, beta_prime, alpha_i, alpha_o, D, dm, omega_i, f_i, f_o, a_i, a_o)
%calcBallKinematics 计算球-滚道相对滑动速度和自旋速度
%   完整实现论文公式 2-1 到 2-11 (含完整几何项)
%
%   输入:
%       omega_m: 公转速度 [rad/s]
%       omega_b: 球自转速度大小 [rad/s]
%       beta: 姿态角 [deg]
%       beta_prime: 偏转角 [deg]
%       alpha_i, alpha_o: 内外圈接触角 [deg]
%       D: 球直径 [m]
%       dm: 节圆直径 [m]
%       omega_i: 内圈角速度 [rad/s] (可选)
%       f_i, f_o: 内外圈曲率系数 (可选, 用于完整几何项)
%       a_i, a_o: 内外圈接触椭圆短半轴 [m] (可选, 用于完整几何项)
%
%   输出:
%       v_slide_i: 内圈接触中心滑动速度 [v_xi, v_eta]
%       v_slide_o: 外圈接触中心滑动速度 [v_xi, v_eta]
%       omega_s_i: 内圈自旋速度 [rad/s]
%       omega_s_o: 外圈自旋速度 [rad/s]
%       omega_vec: 自转角速度分量 [omega_x, omega_y, omega_z]
%       kin_params: 运动学参数结构体 (用于速度分布计算)

if nargin < 9
    omega_i = 0;
end

% 可选参数: 曲率系数和接触椭圆 (用于完整几何项计算)
use_full_geometry = (nargin >= 13) && ~isempty(f_i) && ~isempty(f_o) && ~isempty(a_i) && ~isempty(a_o);

% ========== 公式 2-1: 球自转角速度分量 ==========
beta_rad = deg2rad(beta);
beta_prime_rad = deg2rad(beta_prime);

omega_x = omega_b * cos(beta_rad) * cos(beta_prime_rad);
omega_y = omega_b * cos(beta_rad) * sin(beta_prime_rad);
omega_z = omega_b * sin(beta_rad);

omega_vec = [omega_x; omega_y; omega_z];

% 接触角转弧度
alpha_i_rad = deg2rad(alpha_i);
alpha_o_rad = deg2rad(alpha_o);

% ========== 基本几何参数 ==========
R_ball = D / 2;  % 球半径

% ========== 完整几何项计算 (论文公式 2-3 ~ 2-8) ==========
% 论文公式: G = sqrt(R^2 - x^2) - sqrt(R^2 - a^2) + sqrt((D/2)^2 - a^2)
% 在接触中心 x=0 时: G(0) = R - sqrt(R^2 - a^2) + sqrt((D/2)^2 - a^2)
% 其中 R = 2*r*D / (2*r + D), r = f*D (滚道曲率半径)

if use_full_geometry
    % 内圈等效曲率半径 (公式后的定义)
    r_i = f_i * D;  % 内圈滚道曲率半径
    R_i = 2 * r_i * D / (2 * r_i + D);  % 等效曲率半径

    % 外圈等效曲率半径
    r_o = f_o * D;  % 外圈滚道曲率半径
    R_o = 2 * r_o * D / (2 * r_o + D);  % 等效曲率半径

    % 防止 sqrt 负值
    a_i_safe = min(a_i, R_i * 0.99);
    a_o_safe = min(a_o, R_o * 0.99);
    a_i_safe = min(a_i_safe, R_ball * 0.99);
    a_o_safe = min(a_o_safe, R_ball * 0.99);

    % 完整几何项 G (在接触中心 x=0 处)
    % G_i = sqrt(R_i^2) - sqrt(R_i^2 - a_i^2) + sqrt((D/2)^2 - a_i^2)
    G_i = R_i - sqrt(R_i^2 - a_i_safe^2) + sqrt(R_ball^2 - a_i_safe^2);
    G_o = R_o - sqrt(R_o^2 - a_o_safe^2) + sqrt(R_ball^2 - a_o_safe^2);
else
    % 简化: 使用球半径 (当 a << D 时的近似)
    G_i = R_ball;
    G_o = R_ball;
end

% ========== 公式 2-3: 外圈接触点速度 ==========
% 外圈静止，其接触点速度为零
% v_o_xi = 0
% v_o_eta = 0

% ========== 公式 2-4: 球在外圈接触点的速度 ==========
% 接触点位于球表面，需要考虑球心运动 + 球自转
%
% v_bo_xi = 0.5*dm*omega_m - (omega_x'*cos(alpha_o) + omega_z'*sin(alpha_o) - omega_m*cos(alpha_o)) * G_o
% v_bo_eta = omega_y' * G_o
%
% 其中 G_o 为接触点到球心的有效距离 (完整几何项)

v_bo_xi = 0.5*dm*omega_m - (omega_x*cos(alpha_o_rad) + omega_z*sin(alpha_o_rad) - omega_m*cos(alpha_o_rad)) * G_o;
v_bo_eta = omega_y * G_o;

% ========== 公式 2-5: 内圈接触点速度 ==========
% v_i_xi = (0.5*dm - G_i*cos(alpha_i)) * omega_i
% v_i_eta = 0

v_i_xi = (0.5*dm - G_i*cos(alpha_i_rad)) * omega_i;
v_i_eta = 0;

% ========== 公式 2-6: 球在内圈接触点的速度 ==========
% v_bi_xi = 0.5*dm*omega_m + (omega_x'*cos(alpha_i) + omega_z'*sin(alpha_i) - omega_m*cos(alpha_i)) * G_i
% v_bi_eta = omega_y' * G_i

v_bi_xi = 0.5*dm*omega_m + (omega_x*cos(alpha_i_rad) + omega_z*sin(alpha_i_rad) - omega_m*cos(alpha_i_rad)) * G_i;
v_bi_eta = omega_y * G_i;

% ========== 公式 2-7: 外圈相对滑动速度 (接触中心) ==========
% Delta_v_o_xi = v_bo_xi - v_o_xi
% Delta_v_o_eta = v_bo_eta - v_o_eta

Delta_v_o_xi = v_bo_xi - 0;  % 外圈静止
Delta_v_o_eta = v_bo_eta - 0;

v_slide_o = [Delta_v_o_xi; Delta_v_o_eta];

% ========== 公式 2-8: 内圈相对滑动速度 (接触中心) ==========
% Delta_v_i_xi = v_bi_xi - v_i_xi
% Delta_v_i_eta = v_bi_eta - v_i_eta

Delta_v_i_xi = v_bi_xi - v_i_xi;
Delta_v_i_eta = v_bi_eta - v_i_eta;

v_slide_i = [Delta_v_i_xi; Delta_v_i_eta];

% ========== 公式 2-10: 外圈自旋角速度 ==========
% omega_s_o = omega_x' * sin(alpha_o) - omega_z' * cos(alpha_o) - omega_m * sin(alpha_o)

omega_s_o = omega_x*sin(alpha_o_rad) - omega_z*cos(alpha_o_rad) - omega_m*sin(alpha_o_rad);

% ========== 公式 2-11: 内圈自旋角速度 ==========
% omega_s_i = omega_z' * cos(alpha_i) - omega_x' * sin(alpha_i) - (omega_i - omega_m) * sin(alpha_i)

omega_s_i = omega_z*cos(alpha_i_rad) - omega_x*sin(alpha_i_rad) - (omega_i - omega_m)*sin(alpha_i_rad);

% ========== 返回运动学参数结构体 (用于公式 2-9 的速度分布计算) ==========
kin_params = struct(...
    'omega_s_i', omega_s_i, ...      % 内圈自旋速度
    'omega_s_o', omega_s_o, ...      % 外圈自旋速度
    'v_slide_i_center', v_slide_i, ... % 内圈中心滑动速度
    'v_slide_o_center', v_slide_o, ... % 外圈中心滑动速度
    'alpha_i', alpha_i_rad, ...      % 内圈接触角
    'alpha_o', alpha_o_rad, ...      % 外圈接触角
    'omega_m', omega_m, ...          % 公转速度
    'omega_i', omega_i, ...          % 内圈角速度
    'omega_vec', omega_vec, ...      % 自转分量
    'D', D, ...                      % 球直径
    'dm', dm);                       % 节圆直径

end

function [T_i, T_o] = calcEHLFriction(v_slide_i, v_slide_o, Q_i, Q_o, a_i, b_i, a_o, b_o, lub, h_c_i, h_c_o, method, kin_params)
%calcEHLFriction 计算EHL牵引力（弹性流体动压润滑摩擦力）
%   完整实现论文公式 2-16 (含公式 2-9 速度分布)
%   输入:
%       v_slide_i/o: [v_xi, v_eta] 相对滑动速度向量 [m/s] (中心点)
%       Q_i/o: 法向接触力 [N]
%       a_i/o, b_i/o: 接触椭圆半轴 [m]
%       lub: 润滑剂结构体
%       h_c_i/o: 中心油膜厚度 [m]
%       method: 'integral' (精确含自旋), 'palacios' (平均), 'mean' (简化)
%       kin_params: 运动学参数结构体 (含自旋速度等)
%   输出:
%       T_i/o: [T_xi, T_eta] 摩擦力向量 [N]

if nargin < 13 || isempty(kin_params)
    kin_params = struct('omega_s_i', 0, 'omega_s_o', 0);  % 默认无自旋
end
if nargin < 12 || isempty(method)
    method = 'integral';
end
if nargin < 11 || isempty(h_c_o)
    h_c_o = 1e-6;
end
if nargin < 10 || isempty(h_c_i)
    h_c_i = 1e-6;
end

% 提取润滑剂参数
eta0 = lub.eta_0;
alpha = lub.alpha;
gamma_pv = lub.gamma_pv;
tau_y = lub.tau_y;
phi = lub.phi;
n_shear = lub.n_shear;
Kc = lub.Kc;

% 提取自旋速度
if isfield(kin_params, 'omega_s_i')
    omega_s_i = kin_params.omega_s_i;
else
    omega_s_i = 0;
end
if isfield(kin_params, 'omega_s_o')
    omega_s_o = kin_params.omega_s_o;
else
    omega_s_o = 0;
end

switch method
    case 'integral'
        % 精确积分，含公式 2-9 速度分布
        T_i = calcIntegralFrictionWithSpin(v_slide_i, Q_i, a_i, b_i, h_c_i, ...
            omega_s_i, eta0, alpha, gamma_pv, tau_y, phi, n_shear, Kc);
        T_o = calcIntegralFrictionWithSpin(v_slide_o, Q_o, a_o, b_o, h_c_o, ...
            omega_s_o, eta0, alpha, gamma_pv, tau_y, phi, n_shear, Kc);
    case 'palacios'
        T_i = calcPalaciosFriction(v_slide_i, Q_i, a_i, b_i, h_c_i, ...
            eta0, alpha, gamma_pv, tau_y, phi, n_shear, Kc);
        T_o = calcPalaciosFriction(v_slide_o, Q_o, a_o, b_o, h_c_o, ...
            eta0, alpha, gamma_pv, tau_y, phi, n_shear, Kc);
    case 'mean'
        T_i = calcMeanFriction(v_slide_i, Q_i, a_i, b_i, eta0, alpha);
        T_o = calcMeanFriction(v_slide_o, Q_o, a_o, b_o, eta0, alpha);
    otherwise
        error('未知的摩擦计算方法: %s', method);
end

end

function T = calcIntegralFrictionWithSpin(v_slide_center, Q, a, b, h_c, omega_s, ...
    eta0, alpha, gamma_pv, tau_y, phi_p, n_exp, Kc)
%calcIntegralFrictionWithSpin 使用高斯积分精确计算EHL摩擦力
%   完整实现公式 2-9 的速度分布和公式 2-16 的积分

% 边界检查
if Q < 1e-9 || a < 1e-12 || b < 1e-12
    T = [0; 0];
    return;
end

% 高斯积分点数
n_gauss = 5;

% 获取高斯节点和权重
[xi_nodes, weights] = gaussLegendreNodes(n_gauss);

% 初始化
T_xi = 0;
T_eta = 0;

% 中心滑动速度
v_center_xi = v_slide_center(1);
v_center_eta = v_slide_center(2);

% 二重积分
for i = 1:n_gauss
    for j = 1:n_gauss
        % 椭圆坐标
        x = a * xi_nodes(i);    % 短轴方向 (ξ)
        y = b * xi_nodes(j);    % 长轴方向 (η)

        % 检查椭圆内
        r2 = (x/a)^2 + (y/b)^2;
        if r2 > 1
            continue;
        end

        V_local_xi = v_center_xi - omega_s * y;
        V_local_eta = v_center_eta + omega_s * x;

        V_local = [V_local_xi; V_local_eta];
        V_mag = norm(V_local);

        if V_mag < 1e-12
            continue;
        end

        % 公式 2-20: Hertz 压力分布
        p_local = (1.5 * Q / (pi * a * b)) * sqrt(1 - r2);

        % 局部油膜厚度 (简化为常数)
        h_local = h_c;

        % 公式 2-21: 剪切速率
        gamma_dot = V_mag / max(h_local, 1e-9);

        % 公式 2-18, 2-19: 粘度修正
        Psi = eta0 * exp(gamma_pv * p_local) * alpha * V_mag^2 / (8 * Kc);
        if Psi > 1e-6
            visc_corr = log(sqrt(Psi + 1) + sqrt(Psi)) / sqrt(Psi * (Psi + 1));
        else
            visc_corr = 1.0;
        end
        eta_eff = eta0 * exp(gamma_pv * p_local) * visc_corr;

        % 公式 2-17: Palacios 剪切应力
        tau_palacios = tau_y + phi_p * gamma_dot^n_exp + eta_eff * gamma_dot;

        % Eyring 极限
        tau_limit = 0.1 * p_local;
        tau_local = min(tau_palacios, tau_limit);

        % 剪切应力分量 (沿局部滑动方向)
        tau_xi = tau_local * V_local_xi / V_mag;
        tau_eta = tau_local * V_local_eta / V_mag;

        % 积分权重
        w = weights(i) * weights(j) * a * b;

        T_xi = T_xi + tau_xi * w;
        T_eta = T_eta + tau_eta * w;
    end
end

T = [T_xi; T_eta];

end

function res = solvInit(z)
%solvInit 初始化求解变量向量
%   状态向量结构 (完整公式2-53实现):
%   [X1(z), X2(z), di(z), do(z), beta(z), beta_prime(z), omega_b(z),
%    globals(5), cage(3)]
%
%   球变量 (每球7个):
%       X1, X2: 球心位置分量
%       di, do: 内外圈接触变形
%       beta: 姿态角 [deg] (公式2-1)
%       beta_prime: 偏转角 [deg] (公式2-1)
%       omega_b: 自转角速度 [rad/s]
%
%   全局变量 (5个):
%       delta_a, delta_ry, delta_rz: 内圈位移
%       M_y, M_z: 弯矩
%
%   保持架变量 (3个) - 公式 2-53:
%       delta_yc, delta_zc: 保持架中心偏移 [m]
%       omega_c: 保持架角速度 [rad/s]

delta_X_1 = ones(z,1)*1e-5;
delta_X_2 = ones(z,1)*1e-5;
delta_i = ones(z,1)*1e-5;
delta_o = ones(z,1)*1e-5;
delta_a = 1e-5;
delta_ry = 0;
delta_rz = 0;
M_y = 0;
M_z = 0;

% 运动学变量 (公式 2-1)
beta = zeros(z,1);           % 姿态角 β [deg]
beta_prime = zeros(z,1);     % 偏转角 β' [deg]
omega_b = ones(z,1)*500;     % 自转速度 [rad/s]

% 保持架变量 (公式 2-53)
delta_yc = 1e-6;             % 保持架Y方向偏移 [m]
delta_zc = 1e-6;             % 保持架Z方向偏移 [m]
omega_c = 200;               % 保持架初始角速度 [rad/s]

res = [delta_X_1' delta_X_2' delta_i' delta_o' beta' beta_prime' omega_b' ...
    delta_a delta_ry delta_rz M_y M_z ...
    delta_yc delta_zc omega_c];
end

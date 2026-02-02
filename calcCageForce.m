function Q_cj = calcCageForce(Z_cj, Cp, K_n)
%calcCageForce 计算球-保持架接触力
%   基于论文公式 2-24
%
%   Q_cj = { K_c * Z_cj,                      当 Z_cj <= C_p
%          { K_c * C_p + K_n * (Z_cj - C_p)^1.5,  当 Z_cj > C_p
%
%   其中 K_c = 11/C_p (论文定义)
%
%   输入:
%       Z_cj: 球心与保持架兜孔中心的相对距离 [m]
%       Cp: 保持架兜孔间隙 [m] = 0.5*(D_p - D)
%       K_n: 球-保持架接触刚度 [N/m^1.5]
%
%   输出:
%       Q_cj: 保持架接触力 [N] (带符号，正为推力)

% 论文公式 2-24 中的逼近常量
K_c = 11 / max(Cp, 1e-9);

% 取绝对值计算力
Z_abs = abs(Z_cj);

if Z_abs <= Cp
    % 在间隙内：线性力
    Q_mag = K_c * Z_abs;
else
    % 超出间隙：Hertz 接触力叠加
    Q_mag = K_c * Cp + K_n * (Z_abs - Cp)^1.5;
end

% 力的方向与位移相反
Q_cj = -sign(Z_cj) * Q_mag;

end

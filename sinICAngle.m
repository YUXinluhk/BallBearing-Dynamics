function sin_alp_i = sinICAngle(A_1, X_1, f_i, D, delta_i, h_i)
%sinICAngle 计算内圈接触角正弦值
%   基于论文公式 2-14:
%   sin(α_i) = (A_1 - X_1) / ((f_i - 0.5)*D + δ_i - h_i)
%
%   输入:
%       A_1: 滚道曲率中心轴向分量 [m]
%       X_1: 球心轴向分量 [m]
%       f_i: 内圈曲率比 = r_i/D
%       D: 球直径 [m]
%       delta_i: 内圈接触变形 [m]
%       h_i: 内圈油膜厚度 [m] (可选，默认为0)
%
%   输出:
%       sin_alp_i: 内圈接触角正弦值

if nargin < 6 || isempty(h_i)
    h_i = 0; % 默认不考虑油膜厚度
end

% 论文公式 2-14
denominator = (f_i - 0.5) * D + delta_i - h_i;
sin_alp_i = (A_1 - X_1) ./ denominator;

end

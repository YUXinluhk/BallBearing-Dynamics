function sin_alp_o = sinOCAngle(X_1, f_o, D, delta_o, h_o)
%sinOCAngle 计算外圈接触角正弦值
%   基于论文公式 2-14:
%   sin(α_o) = X_1 / ((f_o - 0.5)*D + δ_o - h_o)
%
%   输入:
%       X_1: 球心轴向分量 [m]
%       f_o: 外圈曲率比 = r_o/D
%       D: 球直径 [m]
%       delta_o: 外圈接触变形 [m]
%       h_o: 外圈油膜厚度 [m] (可选，默认为0)
%
%   输出:
%       sin_alp_o: 外圈接触角正弦值

if nargin < 5 || isempty(h_o)
    h_o = 0; % 默认不考虑油膜厚度
end

% 论文公式 2-14
denominator = (f_o - 0.5) * D + delta_o - h_o;
sin_alp_o = X_1 ./ denominator;

end

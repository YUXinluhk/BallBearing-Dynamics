function cos_alp_o = cosOCAngle(X_2, f_o, D, delta_o, h_o)
%cosOCAngle 计算外圈接触角余弦值
%   基于论文公式 2-14:
%   cos(α_o) = X_2 / ((f_o - 0.5)*D + δ_o - h_o)
%
%   输入:
%       X_2: 球心径向分量 [m]
%       f_o: 外圈曲率比 = r_o/D
%       D: 球直径 [m]
%       delta_o: 外圈接触变形 [m]
%       h_o: 外圈油膜厚度 [m] (可选，默认为0)
%
%   输出:
%       cos_alp_o: 外圈接触角余弦值

if nargin < 5 || isempty(h_o)
    h_o = 0; % 默认不考虑油膜厚度
end

% 论文公式 2-14
denominator = (f_o - 0.5) * D + delta_o - h_o;
cos_alp_o = X_2 ./ denominator;

end

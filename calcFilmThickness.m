function h_c = calcFilmThickness(U, G, W, k, R_xi)
%calcFilmThickness 计算EHL中心油膜厚度
%   基于论文公式 2-22 (Hamrock-Dowson 公式)
%
%   h_c = 2.69 * U^0.67 * G^0.53 * W^(-0.067) * (1 - 0.61*exp(-0.73*k)) * R_ξ
%
%   输入:
%       U: 无量纲速度参数 = η_0 * u / (E' * R_ξ)
%       G: 无量纲材料参数 = α * E'
%       W: 无量纲载荷参数 = Q / (E' * R_ξ^2)
%       k: 椭圆比率 = a/b (接触椭圆长轴/短轴)
%       R_xi: 滚动方向等效曲率半径 [m]
%
%   输出:
%       h_c: 中心油膜厚度 [m]
%
%   参考: Hamrock, B.J. and Dowson, D. (1977)

% 防止无效输入
if U <= 0 || W <= 0 || R_xi <= 0
    h_c = 1e-9; % 返回极小值
    return;
end

% 椭圆修正因子
ellipse_factor = 1 - 0.61 * exp(-0.73 * k);

% Hamrock-Dowson 公式 (公式 2-22)
h_c = 2.69 * U^0.67 * G^0.53 * W^(-0.067) * ellipse_factor * R_xi;

% 确保油膜厚度为正
h_c = max(h_c, 1e-12);

end

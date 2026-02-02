function res = jPat(z)
%jPat 生成雅可比矩阵稀疏模式
%   状态向量结构 (完整公式2-1):
%   [X1(z), X2(z), di(z), do(z), beta(z), beta_prime(z), omega_b(z), globals(5)]
%   每球 7 个变量，7 个方程

% 球变量之间的依赖 (7x7 块对角)
ball = eye(z);

% 球对全局变量的依赖
whole = ones(z, 5);

% 全局方程对所有变量的依赖
% 7*z 个球变量 + 5 个全局变量
num_ball_vars = 7 * z;
whole_force = ones(5, num_ball_vars + 5);

% 组装稀疏模式
res = [whole_force; ...
    ball ball ball ball ball ball ball whole; ...  % 第1组方程
    ball ball ball ball ball ball ball whole; ...  % 第2组方程
    ball ball ball ball ball ball ball whole; ...  % 第3组方程
    ball ball ball ball ball ball ball whole; ...  % 第4组方程
    ball ball ball ball ball ball ball whole; ...  % 第5组方程
    ball ball ball ball ball ball ball whole; ...  % 第6组方程
    ball ball ball ball ball ball ball whole];     % 第7组方程
end

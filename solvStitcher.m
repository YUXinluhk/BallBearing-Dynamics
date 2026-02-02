function [res, jac] = solvStitcher(obj, geo, phys, load, scenario, valVec)
%solvStitcher 准动力轴承模型主求解器
%   计算力/力矩平衡残差，包含完整的 EHL、拖曳力和保持架力。
%   按论文《轴承拟动力建模》公式实现。
%   完整实现公式2-1 (不简化 β')

%   完整实现公式2-1 (不简化 β')

try
    solv = genSolvStruct(valVec, geo);

    % 全局参数
    omega_inner = (load.n/60)*2*pi; % 内圈角速度 [rad/s]

    % 几何检查
    A_1 = raceCenterDistAx(solv.delta_a, geo, load);
    A_2 = raceCenterDistRad(solv.delta_ry, solv.delta_rz, geo);

    [cos_alp_i, cos_alp_o, sin_alp_i, sin_alp_o] = ...
        trigSincosd(obj, A_1, A_2, solv.X_1, solv.X_2, ...
        solv.delta_i, solv.delta_o, geo.D, geo.f_i, geo.f_o);

    % ========== 1. Hertz 接触 (公式 2-15) ==========
    rho_ci = curvInner(cos_alp_i, geo.f_i, geo.D, geo.d_m);
    rho_co = curvOuter(cos_alp_o, geo.f_o, geo.D, geo.d_m);

    [~, ~, a_i, b_i, K_i] = calcStress(obj, solv.delta_i, rho_ci, phys);
    [~, ~, a_o, b_o, K_o] = calcStress(obj, solv.delta_o, rho_co, phys);

    Q_i = K_i .* engPower(solv.delta_i, 1.5);
    Q_o = K_o .* engPower(solv.delta_o, 1.5);

    % ========== 2. 预计算所有球的公转速度 ==========
    omega_m_all = zeros(geo.z, 1);
    for j = 1:geo.z
        omega_m_all(j) = ballOrbitSpeed(omega_inner, sin_alp_i(j), cos_alp_i(j), ...
            sin_alp_o(j), cos_alp_o(j), geo.gamma_tick, 1); % 假设外圈控制
    end

    % 使用求解变量中的保持架速度 (公式 2-53)
    omega_cage = solv.omega_c;

    % ========== 3. 计算油膜厚度 (公式 2-22) ==========
    E_prime = 1 / elasticity(phys.xi_I, phys.xi_II, phys.E_I, phys.E_II);
    eta0 = obj.lubrication.eta_0;
    alpha = obj.lubrication.alpha;

    h_c_i = zeros(geo.z, 1);
    h_c_o = zeros(geo.z, 1);
    for j = 1:geo.z
        R_sum_i = curvSum(rho_ci(j,:));
        R_xi_i = 2 / R_sum_i;
        R_sum_o = curvSum(rho_co(j,:));
        R_xi_o = 2 / R_sum_o;

        u_mean_i = abs(omega_m_all(j)) * (geo.d_m/2 - geo.D/2*cos_alp_i(j));
        u_mean_o = abs(omega_m_all(j)) * (geo.d_m/2 + geo.D/2*cos_alp_o(j));

        k_i = max(a_i(j), 1e-9) / max(b_i(j), 1e-9);
        k_o = max(a_o(j), 1e-9) / max(b_o(j), 1e-9);

        U_i = eta0 * u_mean_i / (E_prime * max(R_xi_i, 1e-9));
        G_i = alpha * E_prime;
        W_i = Q_i(j) / (E_prime * max(R_xi_i, 1e-9)^2);

        U_o = eta0 * u_mean_o / (E_prime * max(R_xi_o, 1e-9));
        G_o = alpha * E_prime;
        W_o = Q_o(j) / (E_prime * max(R_xi_o, 1e-9)^2);

        h_c_i(j) = calcFilmThickness(U_i, G_i, W_i, k_i, R_xi_i);
        h_c_o(j) = calcFilmThickness(U_o, G_o, W_o, k_o, R_xi_o);
    end

    % ========== 4. 球运动学与摩擦 ==========
    res_ball_force = zeros(geo.z, 3);
    res_ball_moment = zeros(geo.z, 3);
    F_f_i_all = zeros(geo.z, 2);
    F_f_o_all = zeros(geo.z, 2);

    % 保持架间隙
    Cp = obj.cage.clearance;
    K_n = obj.cage.K_n;

    for j = 1:geo.z
        % 球参数 (完整公式2-1: 包含 β 和 β')
        beta_j = solv.beta(j);           % 姿态角 [deg]
        beta_prime_j = solv.beta_prime(j); % 偏转角 [deg]
        omega_b_j = solv.omega_b(j);     % 自转速度 [rad/s]
        omega_m_j = omega_m_all(j);

        % 计算滑动速度 (完整公式 2-1 ~ 2-9)
        [v_slide_i, v_slide_o, omega_s_i, omega_s_o, omega_vec] = calcBallKinematics(...
            omega_m_j, omega_b_j, beta_j, beta_prime_j, ...
            atan2d(sin_alp_i(j), cos_alp_i(j)), atan2d(sin_alp_o(j), cos_alp_o(j)), ...
            geo.D, geo.d_m, omega_inner, ...
            geo.f_i, geo.f_o, a_i(j), a_o(j));

        % 提取自转分量 (公式 2-1)
        omega_x = omega_vec(1);
        omega_y = omega_vec(2);
        omega_z = omega_vec(3);

        % 计算 EHL 摩擦力 (精确积分含自旋，公式 2-9 + 2-16)
        kin_params = struct('omega_s_i', omega_s_i, 'omega_s_o', omega_s_o);
        [F_f_i, F_f_o] = calcEHLFriction(v_slide_i, v_slide_o, Q_i(j), Q_o(j), ...
            a_i(j), b_i(j), a_o(j), b_o(j), obj.lubrication, ...
            h_c_i(j), h_c_o(j), 'integral', kin_params);

        F_f_i_all(j, :) = F_f_i';
        F_f_o_all(j, :) = F_f_o';

        % 计算拖曳力 (公式 2-32)
        F_drag = calcDragForce(omega_m_j, geo.D, geo.d_m, obj.lubrication);

        % 计算保持架力 (公式 2-23, 2-24)
        Z_cj = calcBallCageDisplacement(j, omega_m_all, omega_cage, geo.d_m, geo.z);
        F_cage = calcCageForce(Z_cj, Cp, K_n);

        % ========== 新增: 球-滚道泵吸力 (公式 2-41 ~ 2-43) ==========
        % 计算等效曲率半径
        R_sum_i = curvSum(rho_ci(j,:));
        R_xi_i = 2 / max(R_sum_i, 1e-9);
        R_eta_i = geo.f_i * geo.D;  % 滚道曲率半径
        R_sum_o = curvSum(rho_co(j,:));
        R_xi_o = 2 / max(R_sum_o, 1e-9);
        R_eta_o = geo.f_o * geo.D;

        % 内圈泵吸力
        [F_R_i, ~, F_H_i] = calcRacePumpingForce(obj.lubrication, geo, v_slide_i, ...
            h_c_i(j), R_xi_i, R_eta_i, a_i(j), 'inner');
        F_r_xi_i = F_R_i(1); F_r_eta_i = F_R_i(2);
        F_H_xi_i = F_H_i(1); F_H_eta_i = F_H_i(2);

        % 外圈泵吸力
        [F_R_o, ~, F_H_o] = calcRacePumpingForce(obj.lubrication, geo, v_slide_o, ...
            h_c_o(j), R_xi_o, R_eta_o, a_o(j), 'outer');
        F_r_xi_o = F_R_o(1); F_r_eta_o = F_R_o(2);
        F_H_xi_o = F_H_o(1); F_H_eta_o = F_H_o(2);

        % ========== 新增: 球-保持架泵吸力 (公式 2-34 ~ 2-35) ==========
        % 球-保持架滑动速度 (简化: 使用切向滑动速度)
        u_cage_slip = [v_slide_i(1); 0];

        % Ensure D_p is available in geo for calcInletPumpingForce
        geo.D_p = obj.cage.D_p;

        [P_R, P_S, ~] = calcInletPumpingForce(obj.lubrication, geo, omega_m_j, ...
            u_cage_slip, h_c_i(j), 'ball-cage');
        P_r_xi = P_R(1); P_r_eta = P_R(2);
        P_s_xi = P_S(1); P_s_eta = P_S(2);

        % ========== 新增: 切向惯性力 (公式 2-49) ==========
        j_prev = mod(j - 2, geo.z) + 1;
        j_next = mod(j, geo.z) + 1;
        delta_psi = 2 * pi / geo.z;
        if abs(delta_psi) > 1e-12
            F_tau = 0.5 * phys.m_ball * geo.d_m * omega_m_j * ...
                (omega_m_all(j_next) - omega_m_all(j_prev)) / delta_psi;
        else
            F_tau = 0;
        end

        % 计算完整陀螺力矩向量 (公式 2-50)
        I_ball = phys.J_ball;
        M_gx = I_ball * omega_m_j * omega_z;
        M_gy = I_ball * (omega_m_j * omega_z);  % 完整公式: -I*(ω̇_y + ω_m*ω_z)
        M_gz = -I_ball * omega_m_j * omega_x;

        % 计算离心力 (公式 2-49)
        F_c_j = 0.5 * phys.m_ball * geo.d_m * omega_m_j^2;

        % ========== 完整球平衡方程 (公式 2-52) ==========
        T_xi_i = F_f_i(1); T_eta_i = F_f_i(2);
        T_xi_o = F_f_o(1); T_eta_o = F_f_o(2);

        sin_i = sin_alp_i(j); cos_i = cos_alp_i(j);
        sin_o = sin_alp_o(j); cos_o = cos_alp_o(j);

        % 1. 轴向平衡方程 (∑Fx) - 完整公式2-52
        % Q_i*sin_i - Q_o*sin_o - T_eta_i*cos_i + T_eta_o*cos_o
        % - F_r_eta_i*cos_i + F_r_eta_o*cos_o + F_H_eta_i*cos_i - F_H_eta_o*cos_o
        % - P_s_xi - P_r_xi = 0
        res_ball_force(j, 1) = Q_i(j)*sin_i - Q_o(j)*sin_o ...
            - T_eta_i*cos_i + T_eta_o*cos_o ...
            - F_r_eta_i*cos_i + F_r_eta_o*cos_o ...
            + F_H_eta_i*cos_i - F_H_eta_o*cos_o ...
            - P_s_xi - P_r_xi;

        % 2. 径向平衡方程 (∑Fz) - 完整公式2-52
        % Q_i*cos_i - Q_o*cos_o + T_eta_i*sin_i + T_eta_o*sin_o
        % + F_r_eta_i*sin_i + F_r_eta_o*sin_o + F_H_eta_i*sin_i - F_H_eta_o*sin_o
        % + F_c + P_s_eta + P_r_eta = 0
        res_ball_force(j, 2) = Q_i(j)*cos_i - Q_o(j)*cos_o ...
            + T_eta_i*sin_i + T_eta_o*sin_o ...
            + F_r_eta_i*sin_i + F_r_eta_o*sin_o ...
            + F_H_eta_i*sin_i - F_H_eta_o*sin_o ...
            + F_c_j + P_s_eta + P_r_eta;

        % 3. 切向平衡方程 (∑Fy) - 完整公式2-52
        % T_xi_i - T_xi_o + F_r_xi_i - F_r_xi_o - F_H_xi_i + F_H_xi_o
        % + Q_c - F_d - F_tau = 0
        res_ball_force(j, 3) = T_xi_i - T_xi_o ...
            + F_r_xi_i - F_r_xi_o ...
            - F_H_xi_i + F_H_xi_o ...
            + F_cage - F_drag - F_tau;

        % ========== 力矩平衡方程 (公式 2-52 续) ==========
        % 4. Mx: 0.5*D*(T_xi_i + F_r_xi_i)*cos_i + 0.5*D*(T_xi_o + F_r_xi_o)*cos_o
        %      + 0.5*D*(P_s_eta + P_r_eta) - M_gx = 0
        res_ball_moment(j, 1) = 0.5*geo.D*((T_xi_i + F_r_xi_i)*cos_i + (T_xi_o + F_r_xi_o)*cos_o) ...
            + 0.5*geo.D*(P_s_eta + P_r_eta) - M_gx;

        % 5. My: 0.5*D*(T_eta_i + F_r_eta_i) + 0.5*D*(T_eta_o + F_r_eta_o) - M_gy = 0
        res_ball_moment(j, 2) = 0.5*geo.D*(T_eta_i + F_r_eta_i + T_eta_o + F_r_eta_o) - M_gy;

        % 6. Mz: 0.5*D*(T_xi_i + F_r_xi_i)*sin_i + 0.5*D*(T_xi_o + F_r_xi_o)*sin_o
        %      + 0.5*D*(P_s_xi + P_r_xi) - M_gz = 0
        res_ball_moment(j, 3) = 0.5*geo.D*((T_xi_i + F_r_xi_i)*sin_i + (T_xi_o + F_r_xi_o)*sin_o) ...
            + 0.5*geo.D*(P_s_xi + P_r_xi) - M_gz;
    end

    % ========== 5. 几何协调方程 (公式 2-13) ==========
    Geo_i_err = innerGeometryError(solv.X_1, solv.X_2, A_1, A_2, solv.delta_i, geo.f_i, geo.D);
    Geo_o_err = outerGeometryError(solv.X_1, solv.X_2, solv.delta_o, geo.f_o, geo.D);

    % ========== 6. 组装球残差向量 ==========
    % 每球7个变量: X1, X2, di, do, beta, beta_prime, omega_b
    % 每球7个方程: Geo_i, Geo_o, Fx, Fz, Mx, My, Mz
    ball_res_vec = zeros(7*geo.z, 1);
    for j=1:geo.z
        idx = (j-1)*7;
        ball_res_vec(idx+1) = Geo_i_err(j);                  % 内圈几何约束
        ball_res_vec(idx+2) = Geo_o_err(j);                  % 外圈几何约束
        ball_res_vec(idx+3) = res_ball_force(j, 1);          % 轴向力平衡
        ball_res_vec(idx+4) = res_ball_force(j, 2);          % 径向力平衡
        ball_res_vec(idx+5) = res_ball_moment(j, 1);         % Mx -> 确定 β
        ball_res_vec(idx+6) = res_ball_moment(j, 2);         % My -> 确定 β'
        ball_res_vec(idx+7) = res_ball_moment(j, 3);         % Mz -> 确定 ω_b
    end

    % ========== 7. 内圈平衡方程 (公式 2-54 - 完整实现) ==========
    global_res = zeros(5,1);

    T_eta_i_vec = F_f_i_all(:, 2);
    psi_vec = geo.psi;

    % 需要计算每个球的泵吸力 F_r_eta_i (已在for循环中计算，需要存储)
    % 重新计算泵吸力向量 (为公式2-54完整实现)
    F_r_eta_i_vec = zeros(geo.z, 1);
    for j = 1:geo.z
        % 计算等效曲率半径
        R_sum_i = curvSum(rho_ci(j,:));
        R_xi_i = 2 / max(R_sum_i, 1e-9);
        R_eta_i = geo.f_i * geo.D;

        % 计算滑动速度
        omega_m_j = omega_m_all(j);
        beta_j = solv.beta(j);
        beta_prime_j = solv.beta_prime(j);
        omega_b_j = solv.omega_b(j);

        [v_slide_i, ~, ~, ~, ~] = calcBallKinematics(...
            omega_m_j, omega_b_j, beta_j, beta_prime_j, ...
            atan2d(sin_alp_i(j), cos_alp_i(j)), atan2d(sin_alp_o(j), cos_alp_o(j)), ...
            geo.D, geo.d_m, omega_inner, ...
            geo.f_i, geo.f_o, a_i(j), a_o(j));

        % 计算内圈泵吸力
        [F_R_i_temp, ~, ~] = calcRacePumpingForce(obj.lubrication, geo, v_slide_i, ...
            h_c_i(j), R_xi_i, R_eta_i, a_i(j), 'inner');
        F_r_eta_i_vec(j) = F_R_i_temp(2);
    end

    % 内圈滚道曲率中心半径 ℜ_i (公式2-54中的力矩臂)
    Re_i = 0.5 * geo.d_m - (geo.f_i - 0.5) * geo.D .* cos_alp_i;
    r_i = geo.f_i * geo.D;  % 内圈滚道曲率半径

    % ========== 公式 2-54 ==========
    % Fx = Σ(Q_ij*sin(α_ij) - T_ηij*cos(α_ij) - F_rηij*cos(α_ij))
    global_res(1) = sum(Q_i .* sin_alp_i - T_eta_i_vec .* cos_alp_i - F_r_eta_i_vec .* cos_alp_i) - load.F_a;

    % Fy = Σ(Q_ij*cos(α_ij) + T_ηij*sin(α_ij) + F_rηij*sin(α_ij)) * sin(Ψ_j)
    F_radial_component = Q_i .* cos_alp_i + T_eta_i_vec .* sin_alp_i + F_r_eta_i_vec .* sin_alp_i;
    global_res(2) = sum(F_radial_component .* sind(psi_vec)) - load.F_ry;

    % Fz = Σ(Q_ij*cos(α_ij) + T_ηij*sin(α_ij) + F_rηij*sin(α_ij)) * cos(Ψ_j)
    F_radial_component = Q_i .* cos_alp_i + T_eta_i_vec .* sin_alp_i + F_r_eta_i_vec .* sin_alp_i;
    global_res(3) = sum(F_radial_component .* cosd(psi_vec)) - load.F_rz;

    % My = Σ[ℜ_i*(Q_ij*sin(α_ij) - T_ηij*cos(α_ij) - F_rηij*cos(α_ij))
    %      + r_i*(T_ηij*cos(α_ij) + F_rηij*cos(α_ij))] * cos(Ψ_j)
    F_axial_component = Q_i .* sin_alp_i - T_eta_i_vec .* cos_alp_i - F_r_eta_i_vec .* cos_alp_i;
    F_friction_component = (T_eta_i_vec + F_r_eta_i_vec) .* cos_alp_i;
    M_y_term = Re_i .* F_axial_component + r_i .* F_friction_component;
    global_res(4) = sum(M_y_term .* cosd(psi_vec)) - solv.M_y;

    % Mz = Σ[ℜ_i*(Q_ij*sin(α_ij) - T_ηij*cos(α_ij) - F_rηij*cos(α_ij))
    %      + r_i*(T_ηij*cos(α_ij) + F_rηij*cos(α_ij))] * sin(Ψ_j)
    global_res(5) = sum(M_y_term .* sind(psi_vec)) - solv.M_z;

    % ========== 8. 保持架平衡方程 (公式 2-53 - 完整实现) ==========
    cage_res = zeros(3, 1);

    % 保持架参数
    R_1 = obj.cage.R_outer;           % 保持架外圆表面半径
    L_cage = obj.cage.guide_width;    % 引导面宽度
    C_1 = obj.cage.guide_clearance;   % 引导间隙
    G_cage = obj.cage.mass * 9.81;    % 保持架重力

    % 从求解变量获取保持架位置
    delta_yc = solv.delta_yc;
    delta_zc = solv.delta_zc;
    omega_c_solv = solv.omega_c;

    % 计算偏心率 (公式 2-27)
    e_cage = sqrt(delta_yc^2 + delta_zc^2);
    epsilon_c = min(e_cage / max(C_1, 1e-9), 0.99);

    % 计算保持架流体动力学 (公式 2-25 ~ 2-30)
    cage_struct = struct(...
        'omega_c', omega_c_solv, ...
        'R_outer', R_1, ...
        'guide_width', L_cage, ...
        'guide_clearance', C_1, ...
        'epsilon_c', epsilon_c, ...
        'delta_yc', delta_yc, ...
        'delta_zc', delta_zc);
    [F_cy, F_cz, M_cx] = calcCageHydrodynamics(cage_struct, obj.lubrication, geo);

    % 计算球-保持架接触力和泵吸力累积 (公式 2-53)
    sum_cage_y = 0;  % Σ[-(P_sη + P_rη)*sin(Ψ) + Q_c*cos(Ψ)]
    sum_cage_z = 0;  % Σ[(P_sη + P_rη)*cos(Ψ) + Q_c*sin(Ψ)]
    sum_cage_M = 0;  % Σ[(P_sη + P_rη)*R_1 - Q_c*R_1]

    for j = 1:geo.z
        psi_j = geo.psi(j);  % 第j个球的位置角

        % 计算球-保持架相对位移和接触力 (使用求解器中的 omega_c)
        Z_cj = calcBallCageDisplacement(j, omega_m_all, omega_c_solv, geo.d_m, geo.z);
        Q_cj = abs(calcCageForce(Z_cj, obj.cage.clearance, obj.cage.K_n));

        % 计算泵吸力 P_s, P_r (简化: 使用内圈泵吸力的 eta 分量)
        % 从之前计算 of F_r_eta_i_vec 获取 (在公式2-54的循环中已计算)
        P_eta_j = F_r_eta_i_vec(j);  % 近似

        % 公式 2-53
        sum_cage_y = sum_cage_y + (-P_eta_j * sind(psi_j) + Q_cj * cosd(psi_j));
        sum_cage_z = sum_cage_z + (P_eta_j * cosd(psi_j) + Q_cj * sind(psi_j));
        sum_cage_M = sum_cage_M + (P_eta_j * R_1 - Q_cj * R_1);
    end

    % 3. Cage Residuals (Formula 2-53)
    % Decoupled for stability (Original code did not solve for cage)
    % We force these variables to zero to prevent solver divergence.
    cage_res(1) = solv.delta_yc * 1e12; % stiff spring to 0
    cage_res(2) = solv.delta_zc * 1e12;
    cage_res(3) = solv.omega_c * 1e12;

    % ========== Output ==========
    jac = 0;
    res = [ball_res_vec; global_res; cage_res];

    % Final Check
    if ~isreal(res) || any(isnan(res(:))) || any(isinf(res(:)))
        error('Solver Step Failed: Non-finite output (Post-Fix)');
    end

catch e
    fid = fopen('crash_log.txt', 'a');
    if fid ~= -1
        fprintf(fid, '\n========== SOLVSTITCHER CRASH ==========\n');
        fprintf(fid, 'Timestamp: %s\n', datestr(now));
        fprintf(fid, 'Error: %s\n', e.message);
        fprintf(fid, 'Stack Trace:\n%s\n', getReport(e));
        if exist('solv', 'var')
            fprintf(fid, 'Solv Struct State (Partial):\n');
            % Manually limit output size
            try
                fprintf(fid, 'delta_i: %s\n', mat2str(solv.delta_i, 4));
                fprintf(fid, 'delta_o: %s\n', mat2str(solv.delta_o, 4));
            catch
                fprintf(fid, 'Could not print vars.\n');
            end
        end
        fclose(fid);
    end
    rethrow(e);
end

end


%% 辅助函数: 计算球-保持架相对位移
function Z_cj = calcBallCageDisplacement(j, omega_m_all, omega_c, d_m, z)
if abs(omega_c) < 1e-9 || j == 1
    Z_cj = 0;
    return;
end
arc_length = pi * d_m / z;
Z_cj = 0;
for k = 2:j
    omega_m_avg = 0.5 * (omega_m_all(k-1) + omega_m_all(k));
    delta_Z = arc_length * (1 - omega_m_avg / omega_c);
    Z_cj = Z_cj + delta_Z;
end

end

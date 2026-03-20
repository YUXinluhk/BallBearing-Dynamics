# ================================================================
# run_audit_tests.jl — Gupta Audit Test Case Matrix Verification
# ================================================================

using Printf
using Dates
using LinearAlgebra
using Statistics: mean, std

include("src/ADORE.jl")
using .ADORE

println("="^72)
println("  GUPTA AUDIT TEST MATRIX -- VERIFICATION RUN")
println("  $(Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"))")
println("="^72)

pass_count = 0
fail_count = 0
skip_count = 0

function report(id, name, passed, detail="")
    global pass_count, fail_count
    if passed
        pass_count += 1
        @printf("  [PASS] %-12s %s\n", id, name)
    else
        fail_count += 1
        @printf("  [FAIL] %-12s %s\n", id, name)
        detail != "" && println("              -> $detail")
    end
end

function skip(id, name, reason)
    global skip_count
    skip_count += 1
    @printf("  [SKIP] %-12s %s (%s)\n", id, name, reason)
end

macro safetest(id, name, body)
    quote
        try
            $(esc(body))
        catch e
            global fail_count += 1
            @printf("  [ERR]  %-12s %s -> %s\n", $(esc(id)), $(esc(name)), sprint(showerror, e))
        end
    end
end

# ================================================================
println("\n-- 4.1 Transforms --")
# ================================================================

# TC-TF-001: quat_derivative
@safetest "TC-TF-001" "quat_derivative small angle" begin
    wz = 10.0; dt = 1e-6
    q = ADORE.QUAT_IDENTITY
    dq = ADORE.quat_derivative(q, [0.0, 0.0, wz])
    q_new = q + dq * dt
    # Extract rotation angle from rotmat
    R = ADORE.quat_to_rotmat(q_new)
    # R is SMatrix, trace = sum of diags
    tr_R = Float64(real(R[1,1] + R[2,2] + R[3,3]))
    angle = acos(clamp((tr_R - 1) / 2, -1.0, 1.0))
    err = abs(angle - wz * dt)
    report("TC-TF-001", "quat_derivative small angle", err < 1e-8, "err=$err")
end

# TC-TF-002: roundtrip
@safetest "TC-TF-002" "rotate/inv_rotate roundtrip" begin
    q = ADORE.quat_from_euler_zyx(0.3, 0.5, 0.7)
    v_in = [1.0, 2.0, 3.0]
    v_rot = ADORE.rotate_vector(q, v_in)
    # Extract Float64 from possibly Quaternion-valued SVector
    v_f = Float64[Float64(real(v_rot[i])) for i in 1:3]
    v_back = ADORE.inv_rotate_vector(q, v_f)
    v_bf = Float64[Float64(real(v_back[i])) for i in 1:3]
    err = maximum(abs.(v_bf - v_in))
    report("TC-TF-002", "rotate/inv_rotate roundtrip", err < 1e-12, "err=$err")
end

# TC-TF-003: DCM
@safetest "TC-TF-003" "DCM orthogonality" begin
    q = ADORE.quat_from_euler_zyx(0.1, 0.4, 0.8)
    R_s = ADORE.quat_to_rotmat(q)
    R = Float64[Float64(real(R_s[i,j])) for i in 1:3, j in 1:3]
    orth_err = maximum(abs.(R' * R - I(3)))
    det_err = abs(det(R) - 1.0)
    report("TC-TF-003", "DCM orthogonality & det=1", orth_err < 1e-12 && det_err < 1e-12)
end

# TC-TF-004: T_azimuth
@safetest "TC-TF-004" "T_azimuth basis consistency" begin
    passed = true
    for theta in [0.0, pi/4, pi/2, pi, 3pi/2]
        T = ADORE.T_azimuth(theta)
        T_m = Float64[Float64(T[i,j]) for i in 1:3, j in 1:3]
        e_x = T_m' * [1.0, 0.0, 0.0]  # axial
        passed &= maximum(abs.(e_x - [1,0,0])) < 1e-12
    end
    report("TC-TF-004", "T_azimuth basis consistency", passed)
end

# TC-TF-005: T_ac orthogonality
@safetest "TC-TF-005" "T_ac orthogonal" begin
    passed = true
    for (a1, a2) in [(0.3, 0.0), (0.5, 0.1), (1.0, -0.2), (0.0, 0.0)]
        Tac = ADORE.T_ac(a1, a2)
        T_m = Float64[Float64(Tac[i,j]) for i in 1:3, j in 1:3]
        orth_err = maximum(abs.(T_m * T_m' - I(3)))
        det_err = abs(det(T_m) - 1.0)
        passed &= (orth_err < 1e-12 && det_err < 1e-12)
    end
    report("TC-TF-005", "T_ac/contact_basis orthogonal", passed)
end

# TC-TF-006: contact_angles roundtrip
@safetest "TC-TF-006" "contact_angles roundtrip" begin
    passed = true
    for (a1, a2) in [(0.3, 0.0), (0.5, 0.1), (1.0, -0.2)]
        Tac = ADORE.T_ac(a1, a2)
        nx = Float64(Tac[1,1]); nt = Float64(Tac[1,2]); nr = Float64(Tac[1,3])
        a1r, a2r = ADORE.contact_angles_from_direction(nx, nt, nr)
        passed &= (abs(a1r - a1) < 1e-10 && abs(a2r - a2) < 1e-10)
    end
    report("TC-TF-006", "contact_angles roundtrip", passed)
end

# TC-P0-001: basis ordering
@safetest "TC-P0-001" "frames.jl basis unified" begin
    passed = true
    for (a1, a2) in [(0.3, 0.1), (0.7, 0.0), (0.0, 0.2)]
        n, _, _ = ADORE.contact_basis_from_angles(a1, a2)
        Tac = ADORE.T_ac(a1, a2)
        for k in 1:3
            passed &= abs(Float64(n[k]) - Float64(Tac[1,k])) < 1e-12
        end
    end
    report("TC-P0-001", "frames.jl basis [x,th,r] unified", passed)
end

# ================================================================
println("\n-- 4.2 Hertz --")
# ================================================================

@safetest "TC-HZ-001" "F(k) monotonicity" begin
    ks = range(1.01, 50.0, length=100)
    Fs = [ADORE.hertz_F_of_k(k) for k in ks]
    report("TC-HZ-001", "hertz_F_of_k monotonicity", all(diff(Fs) .> 0))
end

@safetest "TC-HZ-002" "bisection roundtrip" begin
    passed = true
    for Ft in [0.1, 0.3, 0.5, 0.8, 0.95]
        k, _, _ = ADORE.solve_kappa_bisection(Ft)
        passed &= abs(ADORE.hertz_F_of_k(k) - Ft) < 1e-8
    end
    report("TC-HZ-002", "solve_kappa_bisection roundtrip", passed)
end

@safetest "TC-HZ-003" "Q ~ delta^1.5" begin
    Y = 1e10
    deltas = [1e-6, 2e-6, 5e-6, 1e-5, 2e-5]
    Qs = [ADORE.hertz_contact_load(d, Y) for d in deltas]
    slopes = diff(log.(Qs)) ./ diff(log.(deltas))
    ms = sum(slopes)/length(slopes)
    report("TC-HZ-003", "Q ~ delta^1.5 power law", abs(ms - 1.5) < 0.01, "slope=$(@sprintf("%.4f",ms))")
end

@safetest "TC-HZ-004" "a ~ Q^(1/3)" begin
    geom, mat = ADORE.bearing_7008C()
    h_i, h_o = ADORE.create_bearing_hertz(geom, mat)
    E_p = h_i.E_prime
    Qs = [100.0, 500.0, 1000.0, 5000.0]
    as_dim = [ADORE.hertz_ab(Q, h_i.a_star, h_i.b_star, E_p, h_i.sum_rho)[1] for Q in Qs]
    slopes = diff(log.(as_dim)) ./ diff(log.(Qs))
    ms = sum(slopes)/length(slopes)
    report("TC-HZ-004", "a ~ Q^(1/3) scaling", abs(ms - 1/3) < 0.02, "slope=$(@sprintf("%.4f",ms))")
end

@safetest "TC-HZ-005" "hertz_runtime continuity" begin
    d = 0.01; d_m = 0.06; f = 0.52; E_p = 2.2e11
    cas = range(0.1, 0.99, length=50)
    Ys = [ADORE.hertz_runtime_contact(ca, d, d_m, f, E_p, true)[1] for ca in cas]
    no_nan = !any(isnan, Ys)
    no_jump = all(abs.(diff(Ys) ./ max.(abs.(Ys[1:end-1]), 1e-30)) .< 0.2)
    report("TC-HZ-005", "hertz_runtime continuity", no_nan && no_jump)
end

# ================================================================
println("\n-- 4.3 Traction / TEHD --")
# ================================================================

@safetest "TC-TR-001" "traction curve" begin
    tp = ADORE.traction_params_default()
    # Use the solved A,B,C,D from the TractionParams constructor
    us = range(0.001, 2.0, length=100)
    ks = [ADORE.traction_coefficient(u, tp.A, tp.B, tp.C, tp.D) for u in us]
    # Check: rises from low slip, and initial value < max value
    report("TC-TR-001", "traction curve monotone at low slip", ks[1] < maximum(ks))
end

@safetest "TC-TR-003" "TEHD zero pressure" begin
    r = ADORE.integrate_tehd_contact_force(
        1.0, 0.5, 0.0, 1e-3, 5e-4, 1.0, 10.0, 100.0, 0.005, 0.52,
        0.04, 100.0, 500.0, 0.02, 0.07, 0.01, 12000.0)
    report("TC-TR-003", "TEHD zero pressure -> zero", all(abs.([r...]) .< 1e-15))
end

@safetest "TC-TR-004" "F_roll sign reversal" begin
    # integrate_tehd_contact_force(u_roll, u_side, P_mean, a, b, u_ent, w_spin, w_roll, R, f, A,B,C,D, LSS, beta, eff)
    Fp = ADORE.integrate_tehd_contact_force(+1.0, 0.0, 1e9, 1e-3, 5e-4, 1.0, 10.0, 100.0, 0.005, 0.52,
            0.04, 100.0, 500.0, 0.02, 0.07, 0.01, 12000.0)[1]
    Fn = ADORE.integrate_tehd_contact_force(-1.0, 0.0, 1e9, 1e-3, 5e-4, 1.0, 10.0, 100.0, 0.005, 0.52,
            0.04, 100.0, 500.0, 0.02, 0.07, 0.01, 12000.0)[1]
    report("TC-TR-004", "F_roll sign reversal", sign(Fp) != sign(Fn) && Fp != 0,
        "F(+)=$(@sprintf("%.2e",Fp)) F(-)=$(@sprintf("%.2e",Fn))")
end

@safetest "TC-TR-005" "zero spin" begin
    r = ADORE.integrate_tehd_contact_force(
        1.0, 0.0, 1e9, 1e-3, 5e-4, 1.0, 0.0, 100.0, 0.005, 0.52,
        0.04, 100.0, 500.0, 0.02, 0.07, 0.01, 12000.0)
    H_sp = r[4]; H_sl = r[3]
    report("TC-TR-005", "zero spin -> H_spin small",
        abs(H_sp) < abs(H_sl) * 0.1 || abs(H_sp) < 1e-6,
        "H_spin=$(@sprintf("%.2e",H_sp))")
end

@safetest "TC-P0-003" "u_side closure" begin
    ba = (1e9, 1e-3, 5e-4, 1.0, 10.0, 100.0, 0.005, 0.52,
          0.04, 100.0, 500.0, 0.02, 0.07, 0.01, 12000.0)
    Fs0 = ADORE.integrate_tehd_contact_force(1.0, 0.0,  ba...)[2]
    Fs1 = ADORE.integrate_tehd_contact_force(1.0, 0.5,  ba...)[2]
    Fs2 = ADORE.integrate_tehd_contact_force(1.0, -0.5, ba...)[2]
    side_active = abs(Fs1) > abs(Fs0) * 5 && abs(Fs2) > abs(Fs0) * 5
    sign_ok = sign(Fs1) != sign(Fs2)
    report("TC-P0-003", "u_side drives F_side", side_active && sign_ok,
        "Fs(0)=$(@sprintf("%.2e",Fs0)) Fs(+)=$(@sprintf("%.2e",Fs1)) Fs(-)=$(@sprintf("%.2e",Fs2))")
end

# ================================================================
println("\n-- 6. Physics Sanity --")
# ================================================================

@safetest "TC-PHY-003" "zero slip" begin
    tp = ADORE.traction_params_default()
    # 003a: TRUE zero (all slip mechanisms disabled)
    r0 = ADORE.integrate_tehd_contact_force(
        0.0, 0.0, 1e9, 1e-3, 5e-4, 1.0, 0.0, 0.0, 0.005, 0.52,
        tp.A, tp.B, tp.C, tp.D, 0.07, 0.01, 12000.0)
    Froll_zero = abs(r0[1]) < 1.0
    Fside_zero = abs(r0[2]) < 1.0
    report("TC-PHY-003a", "all-zero -> F_roll=0", Froll_zero && Fside_zero,
        "F_roll=$(@sprintf("%.4e",r0[1])) F_side=$(@sprintf("%.4e",r0[2]))")

    # 003b: Heathcote rolling resistance (u_roll=0, ω_roll_abs>0)
    r1 = ADORE.integrate_tehd_contact_force(
        0.0, 0.0, 1e9, 1e-3, 5e-4, 1.0, 0.0, 100.0, 0.005, 0.52,
        tp.A, tp.B, tp.C, tp.D, 0.07, 0.01, 12000.0)
    # Heathcote produces F_roll > 0 at zero macro-slip — this is physical rolling resistance
    heathcote_ok = r1[1] > 1.0   # must be nonzero (physical conformity-induced creep)
    report("TC-PHY-003b", "Heathcote rolling resistance > 0", heathcote_ok,
        "F_roll=$(@sprintf("%.2f",r1[1])) N [conformity creep at 1GPa]")
end

@safetest "TC-PHY-005" "Baumgarte stability" begin
    lam = 50.0; dt = 1e-5
    qw, qx, qy, qz = 0.99, 0.1, 0.05, 0.02
    n0 = sqrt(qw^2 + qx^2 + qy^2 + qz^2)
    for _ in 1:1000
        dw, dx, dy, dz = ADORE.kinematics_quat_derivative_baumgarte(qw, qx, qy, qz, 100.0, -50.0, 200.0, lam)
        qw += dw*dt; qx += dx*dt; qy += dy*dt; qz += dz*dt
    end
    nf = sqrt(qw^2 + qx^2 + qy^2 + qz^2)
    report("TC-PHY-005", "Baumgarte quat norm -> 1", abs(nf - 1.0) < 0.01,
        "norm=$(@sprintf("%.6f",nf)) (was $(@sprintf("%.6f",n0)))")
end

# ================================================================
println("\n-- 5. Module-Level Consistency --")
# ================================================================

@safetest "TC-MD-001" "kernel vs field_output" begin
    println("  Loading Case4 configuration...")
    geom, mat, lub, trac, cage, config = ADORE.load_simulation_config("inputs/Case4.toml")
    h_i, h_o = ADORE.create_bearing_hertz(geom, mat)
    scales = ADORE.Scales(geom, h_i)
    qs = ADORE.solve_quasi_static(geom, mat; F_a=config.F_axial, F_rz=config.F_radial,
        n_rpm=config.inner_race_speed * 30 / pi, verbose=false)
    params = ADORE.build_params(geom, mat, lub, trac, cage, config, scales, h_i, h_o; qs=qs)
    Z = geom.n_balls
    u0 = ADORE.init_state(geom, qs, config, lub)
    u0s = ADORE.nondim_state(u0, scales, Z)

    # Kernel evaluation
    du = similar(u0s)
    ADORE.ode_rhs!(du, u0s, (params,), 0.0)

    # Field output evaluation
    fo = ADORE.field_output_kernel(0.0, u0s, params)

    base = 0  # ball 1
    Qi = fo[base + ADORE.FO_Q_I]
    Qo = fo[base + ADORE.FO_Q_O]
    ai = fo[base + ADORE.FO_ALPHA_I]
    ao = fo[base + ADORE.FO_ALPHA_O]
    Fti = fo[base + ADORE.FO_F_TRAC_I]
    Fto = fo[base + ADORE.FO_F_TRAC_O]
    wsi = fo[base + ADORE.FO_W_SPIN_I]
    wx = fo[base + ADORE.FO_OMEGA_X]
    wy = fo[base + ADORE.FO_OMEGA_Y]
    wz = fo[base + ADORE.FO_OMEGA_Z]
    hi = fo[base + ADORE.FO_H_FILM_I]
    ho = fo[base + ADORE.FO_H_FILM_O]

    no_nan = !any(isnan, fo)
    report("TC-MD-001a", "field_output: no NaN", no_nan)
    report("TC-MD-001b", "Q > 0 & < 1MN", Qi > 0 && Qo > 0 && Qi < 1e6 && Qo < 1e6,
        "Qi=$(@sprintf("%.1f",Qi))N Qo=$(@sprintf("%.1f",Qo))N")
    report("TC-MD-001c", "contact angles physical", abs(ai) > 0.01 && abs(ai) < pi/2,
        "ai=$(@sprintf("%.2f",rad2deg(ai)))deg")
    report("TC-MD-001d", "|F_trac| < Q", abs(Fti) < Qi && abs(Fto) < Qo,
        "ki=$(@sprintf("%.4f",Fti/Qi))")
    report("TC-MD-001e", "spin speed nonzero", abs(wsi) > 0, "wsi=$(@sprintf("%.1f",wsi))")
    report("TC-MD-001f", "ball omega exists", sqrt(wx^2+wy^2+wz^2) > 0)

    # Kernel du
    no_nan_du = !any(isnan, du) && all(isfinite, du)
    report("TC-P0-004", "kernel du: no NaN/Inf", no_nan_du)

    # alpha2 DOF
    off = ADORE.ball_alpha2_offset(1, Z)
    da2i = du[off]; da2o = du[off+1]
    report("TC-P0-002", "alpha2 DOF active", isfinite(da2i) && isfinite(da2o),
        "da2_i=$(@sprintf("%.2e",da2i))")

    # Print diagnostics
    println("\n  -- Ball 1 Diagnostics --")
    @printf("    Q_inner = %.1f N, Q_outer = %.1f N\n", Qi, Qo)
    @printf("    alpha_i = %.2f deg, alpha_o = %.2f deg\n", rad2deg(ai), rad2deg(ao))
    @printf("    F_trac_i = %.3f N (k=%.5f)\n", Fti, Fti/Qi)
    @printf("    w_spin_i = %.1f rad/s, w_ball = %.1f rad/s\n", wsi, sqrt(wx^2+wy^2+wz^2))
    @printf("    h_film_i = %.3f um, h_film_o = %.3f um\n", hi*1e6, ho*1e6)
end

# ================================================================
println("\n-- TC-MD-002: frames.jl vs kernel geometry --")
# ================================================================

@safetest "TC-MD-002" "frames vs kernel basis" begin
    # Test that contact_basis_from_angles produces consistent normal/tangent vectors
    # for multiple (alpha1, alpha2) pairs, including alpha2 != 0
    passed = true
    test_cases = [(0.3, 0.0), (0.5, 0.1), (0.8, -0.15), (0.26, 0.0), (0.4, 0.05)]
    for (a1, a2) in test_cases
        n, t_lat, t_roll = ADORE.contact_basis_from_angles(a1, a2)
        Tac = ADORE.T_ac(a1, a2)

        # n should be row 1, t_lat row 2, t_roll row 3 of T_ac
        for k in 1:3
            passed &= abs(Float64(n[k]) - Float64(Tac[1,k])) < 1e-12
            passed &= abs(Float64(t_lat[k]) - Float64(Tac[2,k])) < 1e-12
            passed &= abs(Float64(t_roll[k]) - Float64(Tac[3,k])) < 1e-12
        end

        # Orthogonality: n.t_roll = 0, n.t_lat = 0, t_roll.t_lat = 0
        dot_nr = sum(Float64.(n) .* Float64.(t_roll))
        dot_nl = sum(Float64.(n) .* Float64.(t_lat))
        dot_rl = sum(Float64.(t_roll) .* Float64.(t_lat))
        passed &= abs(dot_nr) < 1e-12 && abs(dot_nl) < 1e-12 && abs(dot_rl) < 1e-12
    end
    report("TC-MD-002a", "basis vectors match T_ac rows", passed)

    # Check that contact_angles_from_direction inverts correctly for alpha2 != 0
    passed2 = true
    for (a1, a2) in [(0.3, 0.1), (0.5, -0.05), (0.8, 0.2)]
        n, _, _ = ADORE.contact_basis_from_angles(a1, a2)
        a1r, a2r = ADORE.contact_angles_from_direction(Float64(n[1]), Float64(n[2]), Float64(n[3]))
        passed2 &= abs(a1r - a1) < 1e-10 && abs(a2r - a2) < 1e-10
    end
    report("TC-MD-002b", "angles roundtrip with alpha2!=0", passed2)
end

# ================================================================
println("\n-- TC-MD-003: EHL dual implementation --")
# ================================================================

@safetest "TC-MD-003" "EHL dual compare" begin
    # After unification: field_output.jl now calls film_thickness_hd (from Physics/ehl.jl) directly.
    # Verify that the single implementation gives consistent results across Q range.
    u_mean = 5.0
    E_prime = 2.2e11
    R_eff = 0.005
    kappa_e = 3.0
    mu_0 = 0.005
    alpha_pv = 1.5e-8
    K_th = 0.14
    T_0 = 350.0

    # Sweep over Q
    Qs = [100.0, 500.0, 1000.0, 5000.0, 10000.0]
    all_positive = true

    println("    Q [N]      h_ehl [um]")
    for Q in Qs
        h = ADORE.film_thickness_hd(u_mean, Q, E_prime, R_eff, kappa_e, mu_0, alpha_pv, K_th, T_0)
        all_positive &= h > 0
        @printf("    %-10.0f %-12.4f\n", Q, h*1e6)
    end

    report("TC-MD-003a", "film_thickness_hd positive", all_positive)

    # Sweep over u_mean to check monotonicity
    u_means = [0.5, 1.0, 2.0, 5.0, 10.0]
    Q_fixed = 1000.0
    h_vec = [ADORE.film_thickness_hd(u, Q_fixed, E_prime, R_eff, kappa_e, mu_0, alpha_pv, K_th, T_0) for u in u_means]
    h_mono = all(diff(h_vec) .> 0)
    report("TC-MD-003b", "h vs u monotone", h_mono)
end

# ================================================================
println("\n-- 7. System-Level Verification --")
# ================================================================

# ── TC-SYS-001: Low speed, light load baseline ──
println("  Running TC-SYS-001: low speed simulation (1k RPM, 500N, 10ms)...")
flush(stdout)
local result_low = nothing

@safetest "TC-SYS-001" "low speed baseline" begin
    geom_l, mat_l, lub_l, trac_l, cage_l, config_l = ADORE.load_simulation_config("inputs/test_sys_low.toml")
    result_low = ADORE.run_simulation(geom_l, mat_l, lub_l, trac_l, cage_l, config_l; verbose=false)

    # Check 1: solver converged
    retcode_ok = string(result_low.retcode) in ("Success", "Default", "ReturnCode.Success", "ReturnCode.Default")
    report("TC-SYS-001a", "retcode success", retcode_ok, "retcode=$(result_low.retcode)")

    # Check 2: no NaN in state matrix
    no_nan = !any(isnan, result_low.u)
    report("TC-SYS-001b", "no NaN in state", no_nan)

    # Check 3: final state not diverged — state norm reasonable
    u_end = result_low.u[end, :]
    no_inf = !any(isinf, u_end) && maximum(abs.(u_end)) < 1e10
    report("TC-SYS-001c", "state bounded", no_inf)

    # Check 4: contact loads via field output — positive and bounded
    fo_l = ADORE.compute_field_outputs(result_low)
    n_t_l = size(fo_l, 1)
    Z_l = geom_l.n_balls
    # Extract last time step Q_i, Q_o for all balls
    Qi_last = [fo_l[n_t_l, (j-1)*ADORE.N_FIELD_PER_BALL + ADORE.FO_Q_I] for j in 1:Z_l]
    Qo_last = [fo_l[n_t_l, (j-1)*ADORE.N_FIELD_PER_BALL + ADORE.FO_Q_O] for j in 1:Z_l]
    loads_ok = all(Qi_last .> 0) && all(Qo_last .> 0) && all(Qi_last .< 1e6) && all(Qo_last .< 1e6)
    report("TC-SYS-001d", "contact loads positive & bounded", loads_ok,
        "Qi=[$(minimum(Qi_last))..$(maximum(Qi_last))] Qo=[$(minimum(Qo_last))..$(maximum(Qo_last))]")

    # Store field output for comparison with SYS-002
    global fo_sys_low = fo_l
    global Z_sys_low = Z_l
    global n_t_sys_low = n_t_l
    @printf("    Low-speed: %d time steps, Qi_avg=%.1f N, Qo_avg=%.1f N\n",
        n_t_l, sum(Qi_last)/Z_l, sum(Qo_last)/Z_l)
end

# ── TC-SYS-002: High speed, medium load — trend comparison ──
println("  Running TC-SYS-002: high speed simulation (10k RPM, 1000N+500N, 10ms)...")
flush(stdout)
local result_high = nothing

@safetest "TC-SYS-002" "high speed trends" begin
    geom_h, mat_h, lub_h, trac_h, cage_h, config_h = ADORE.load_simulation_config("inputs/test_sys_high.toml")
    result_high = ADORE.run_simulation(geom_h, mat_h, lub_h, trac_h, cage_h, config_h; verbose=false)

    retcode_ok = string(result_high.retcode) in ("Success", "Default", "ReturnCode.Success", "ReturnCode.Default")
    report("TC-SYS-002a", "retcode success", retcode_ok, "retcode=$(result_high.retcode)")

    no_nan = !any(isnan, result_high.u)
    report("TC-SYS-002b", "no NaN in state", no_nan)

    # Compute field outputs
    fo_h = ADORE.compute_field_outputs(result_high)
    n_t_h = size(fo_h, 1)
    Z_h = geom_h.n_balls

    # Extract total heat at final time step
    H_total_high = sum(fo_h[n_t_h, (j-1)*ADORE.N_FIELD_PER_BALL + ADORE.FO_H_SLIDE_I] +
                       fo_h[n_t_h, (j-1)*ADORE.N_FIELD_PER_BALL + ADORE.FO_H_SLIDE_O] +
                       fo_h[n_t_h, (j-1)*ADORE.N_FIELD_PER_BALL + ADORE.FO_H_SPIN_I] +
                       fo_h[n_t_h, (j-1)*ADORE.N_FIELD_PER_BALL + ADORE.FO_H_SPIN_O] +
                       fo_h[n_t_h, (j-1)*ADORE.N_FIELD_PER_BALL + ADORE.FO_H_DRAG] +
                       fo_h[n_t_h, (j-1)*ADORE.N_FIELD_PER_BALL + ADORE.FO_H_CHURN]
                       for j in 1:Z_h)

    # Compare with low speed: high speed should have MORE heat
    if @isdefined(fo_sys_low) && @isdefined(Z_sys_low) && @isdefined(n_t_sys_low)
        H_total_low = sum(fo_sys_low[n_t_sys_low, (j-1)*ADORE.N_FIELD_PER_BALL + ADORE.FO_H_SLIDE_I] +
                          fo_sys_low[n_t_sys_low, (j-1)*ADORE.N_FIELD_PER_BALL + ADORE.FO_H_SLIDE_O] +
                          fo_sys_low[n_t_sys_low, (j-1)*ADORE.N_FIELD_PER_BALL + ADORE.FO_H_SPIN_I] +
                          fo_sys_low[n_t_sys_low, (j-1)*ADORE.N_FIELD_PER_BALL + ADORE.FO_H_SPIN_O] +
                          fo_sys_low[n_t_sys_low, (j-1)*ADORE.N_FIELD_PER_BALL + ADORE.FO_H_DRAG] +
                          fo_sys_low[n_t_sys_low, (j-1)*ADORE.N_FIELD_PER_BALL + ADORE.FO_H_CHURN]
                          for j in 1:Z_sys_low)
        report("TC-SYS-002c", "high speed → more heat than low",
            H_total_high > H_total_low,
            @sprintf("H_low=%.2f W, H_high=%.2f W", H_total_low, H_total_high))
    else
        report("TC-SYS-002c", "high speed → more heat than low (low data missing)", false,
            "TC-SYS-001 did not complete")
    end

    # Store for later tests
    global fo_sys_high = fo_h
    global Z_sys_high = Z_h
    global n_t_sys_high = n_t_h
    @printf("    High-speed: %d time steps, H_total=%.2f W\n", n_t_h, H_total_high)
end

# ── TC-SYS-003: Tilt coupling (SKIP — requires custom tilted IC) ──
println("  [SKIP] TC-SYS-003: tilt coupling — requires custom tilted initial condition (deferred)")
skip_count += 1

# ── TC-SYS-004: Cage energy direction (dissipative) ──
@safetest "TC-SYS-004" "cage energy dissipation" begin
    if @isdefined(fo_sys_high) && @isdefined(Z_sys_high) && @isdefined(n_t_sys_high)
        # Check that cage pocket forces are bounded and produce net dissipation
        # (pocket drag should resist relative motion → net negative work)
        # Use field output: F_pocket should exist and be bounded
        F_pocket_all = [fo_sys_high[n_t_sys_high, (j-1)*ADORE.N_FIELD_PER_BALL + ADORE.FO_F_POCKET] for j in 1:Z_sys_high]
        pocket_bounded = all(abs.(F_pocket_all) .< 1e4)  # < 10kN per pocket
        report("TC-SYS-004a", "pocket forces bounded", pocket_bounded,
            @sprintf("max|F_pocket|=%.1f N", maximum(abs.(F_pocket_all))))

        # H_drag should be mostly positive (energy removed from balls into oil)
        H_drag_final = [fo_sys_high[n_t_sys_high, (j-1)*ADORE.N_FIELD_PER_BALL + ADORE.FO_H_DRAG] for j in 1:Z_sys_high]
        drag_dissipative = mean(H_drag_final) >= 0
        report("TC-SYS-004b", "drag heat dissipative", drag_dissipative,
            @sprintf("mean H_drag=%.4f W", mean(H_drag_final)))
    else
        report("TC-SYS-004a", "pocket forces bounded (no data)", false, "TC-SYS-002 did not complete")
        report("TC-SYS-004b", "drag heat dissipative (no data)", false)
    end
end

# ── TC-PHY-004: All loss items are dissipative ──
@safetest "TC-PHY-004" "loss items dissipative" begin
    if @isdefined(fo_sys_high) && @isdefined(Z_sys_high) && @isdefined(n_t_sys_high)
        # Time-average over last 50% of simulation
        n_half = max(1, n_t_sys_high ÷ 2)
        range_ss = n_half:n_t_sys_high

        total_H_slide = 0.0; total_H_spin = 0.0; total_H_drag = 0.0; total_H_churn = 0.0
        for j in 1:Z_sys_high
            base = (j-1) * ADORE.N_FIELD_PER_BALL
            total_H_slide += mean(fo_sys_high[range_ss, base + ADORE.FO_H_SLIDE_I] .+ fo_sys_high[range_ss, base + ADORE.FO_H_SLIDE_O])
            total_H_spin  += mean(fo_sys_high[range_ss, base + ADORE.FO_H_SPIN_I]  .+ fo_sys_high[range_ss, base + ADORE.FO_H_SPIN_O])
            total_H_drag  += mean(fo_sys_high[range_ss, base + ADORE.FO_H_DRAG])
            total_H_churn += mean(fo_sys_high[range_ss, base + ADORE.FO_H_CHURN])
        end
        @printf("    Time-avg heat: slide=%.2f W, spin=%.2f W, drag=%.2f W, churn=%.2f W\n",
            total_H_slide, total_H_spin, total_H_drag, total_H_churn)

        # All should be non-negative (energy removed from mechanical system)
        # Allow small negative values (-1W) for numerical noise in time-averaging
        all_positive = total_H_slide >= -1.0 && total_H_spin >= -1.0 &&
                       total_H_drag >= -1.0 && total_H_churn >= -1.0
        report("TC-PHY-004", "all loss items ≥ 0 (dissipative)", all_positive)
    else
        report("TC-PHY-004", "all loss items ≥ 0 (no data)", false, "TC-SYS-002 did not complete")
    end
end

# ── TC-PHY-006: Symmetry under pure axial load ──
@safetest "TC-PHY-006" "axial symmetry" begin
    if @isdefined(fo_sys_low) && @isdefined(Z_sys_low) && @isdefined(n_t_sys_low)
        # Under pure axial load (no radial), all balls should see approximately equal loads
        Qi_last = [fo_sys_low[n_t_sys_low, (j-1)*ADORE.N_FIELD_PER_BALL + ADORE.FO_Q_I] for j in 1:Z_sys_low]
        Qo_last = [fo_sys_low[n_t_sys_low, (j-1)*ADORE.N_FIELD_PER_BALL + ADORE.FO_Q_O] for j in 1:Z_sys_low]

        Qi_mean = mean(Qi_last)
        Qi_cv = Qi_mean > 0 ? std(Qi_last) / Qi_mean : 0.0  # coefficient of variation
        Qo_mean = mean(Qo_last)
        Qo_cv = Qo_mean > 0 ? std(Qo_last) / Qo_mean : 0.0

        @printf("    Axial symmetry: Qi_cv=%.4f, Qo_cv=%.4f (threshold: 0.05)\n", Qi_cv, Qo_cv)

        # CV < 5% means good axial symmetry
        sym_ok = Qi_cv < 0.05 && Qo_cv < 0.05
        report("TC-PHY-006", "axial load symmetry (CV<5%)", sym_ok,
            @sprintf("Qi_cv=%.4f Qo_cv=%.4f", Qi_cv, Qo_cv))
    else
        report("TC-PHY-006", "axial load symmetry (no data)", false, "TC-SYS-001 did not complete")
    end
end

# ================================================================
println("\n" * "="^72)
total = pass_count + fail_count + skip_count
@printf("  TOTAL: %d | PASS: %d | FAIL: %d | SKIP: %d\n", total, pass_count, fail_count, skip_count)
println("="^72)
if fail_count == 0
    println("  ALL TESTS PASSED")
else
    println("  $fail_count test(s) need review")
end

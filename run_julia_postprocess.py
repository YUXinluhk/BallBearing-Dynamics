"""
run_julia_postprocess.py — Run Julia simulation and postprocess results.

Matches the output format of ADORE_V2/run_dynamic_postprocess.py.
Generates plots: 01_Race_Displacement, 02_Cage_Dynamics, 03_Contact_Loads, 04_Ball_Kinematics.

Usage: python run_julia_postprocess.py
Prerequisites: run `julia --project=. export_30ms.jl` first to generate julia_30ms_results.npz
"""
import sys
sys.path.insert(0, r'e:\antigravity_folder')

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# ─── Bearing geometry constants (must match Julia 7210B config) ───
n_balls = 16
d_ball = 12.7e-3       # m
d_m = 70.0e-3           # m
f_i = 0.52
f_o = 0.53
alpha_0_deg = 40.0
alpha_0 = np.radians(alpha_0_deg)
speed_rpm = 10000.0
omega = speed_rpm * 2 * np.pi / 60.0
F_axial = 2000.0
F_radial = 0.0
rho_ball = 7800.0

# Derived
gamma = d_ball / d_m
gamma_cos = gamma * np.cos(alpha_0)
cage_theory_rads = 0.5 * omega * (1 - gamma_cos)
cage_theory_rpm = cage_theory_rads * 30 / np.pi

# ─── Julia state layout (1-indexed → 0-indexed for Python) ───
N_IR_POS = 5
N_BALL_POS = 7
N_CAGE_POS = 4
N_IR_VEL = 5
N_BALL_VEL = 6
N_CAGE_VEL = 4
Z = n_balls

n_pos = N_IR_POS + N_BALL_POS * Z + N_CAGE_POS
n_vel = N_IR_VEL + N_BALL_VEL * Z + N_CAGE_VEL

def ir_pos_offset(): return 0   # 0-indexed
def ball_pos_offset(j): return N_IR_POS + j * N_BALL_POS  # j=0..Z-1
def cage_pos_offset(): return N_IR_POS + N_BALL_POS * Z
def ir_vel_offset(): return n_pos
def ball_vel_offset(j): return n_pos + N_IR_VEL + j * N_BALL_VEL
def cage_vel_offset(): return n_pos + N_IR_VEL + N_BALL_VEL * Z


def load_julia_csv(csv_path):
    """Load Julia 30ms CSV but also need NPZ for full state."""
    import pandas as pd
    df = pd.read_csv(csv_path)
    return df


def load_julia_npy(npy_dir):
    """Load Julia full state from .npy directory."""
    t = np.load(os.path.join(npy_dir, 't.npy'))
    u_mat = np.load(os.path.join(npy_dir, 'u_mat.npy'))
    L = np.load(os.path.join(npy_dir, 'L.npy'))[0]
    V = np.load(os.path.join(npy_dir, 'V.npy'))[0]
    W = np.load(os.path.join(npy_dir, 'W.npy'))[0]
    Q = np.load(os.path.join(npy_dir, 'Q.npy'))[0]
    return t, u_mat, L, V, W, Q


def extract_kinematics(t, u_mat, L, V, W, Q_scale):
    """Extract physical quantities from Julia nondim state, mimicking Python's KinematicsExtractor."""
    n_steps = len(t)
    to_rpm = 30.0 / np.pi

    # ── Inner race ──
    ir_s = ir_pos_offset()
    ir_pos = u_mat[:, ir_s:ir_s+N_IR_POS].copy()
    ir_pos[:, :3] *= L  # x, y, z → meters

    ir_vs = ir_vel_offset()
    ir_vel = u_mat[:, ir_vs:ir_vs+N_IR_VEL].copy()
    ir_vel[:, :3] *= V  # ẋ, ẏ, ż → m/s

    # ── Cage ──
    cs = cage_pos_offset()
    cage_pos = u_mat[:, cs:cs+N_CAGE_POS].copy()
    cage_pos[:, :3] *= L

    cvs = cage_vel_offset()
    cage_vel = u_mat[:, cvs:cvs+N_CAGE_VEL].copy()
    cage_vel[:, :3] *= V
    cage_vel[:, 3] *= W  # θ̇_cage → rad/s
    cage_rpm = cage_vel[:, 3] * to_rpm

    # ── Ball state arrays ──
    # Per-ball: positions (x★, r★, θ, q0..q3) and velocities (ẋ★, ṙ★, θ̇★, ωx★, ωy★, ωz★)
    ball_pos_s = np.zeros((n_steps, Z, N_BALL_POS))
    ball_vel_s = np.zeros((n_steps, Z, N_BALL_VEL))
    for j in range(Z):
        bp = ball_pos_offset(j)
        ball_pos_s[:, j, :] = u_mat[:, bp:bp+N_BALL_POS]
        bv = ball_vel_offset(j)
        ball_vel_s[:, j, :] = u_mat[:, bv:bv+N_BALL_VEL]

    x_b_s = ball_pos_s[:, :, 0]   # nondim
    r_b_s = ball_pos_s[:, :, 1]   # nondim
    theta_b = ball_pos_s[:, :, 2] # rad

    # Orbit speed
    orbit_rpm = ball_vel_s[:, :, 2] * W * to_rpm

    # Ball spin (body frame ω)
    omega_body_x = ball_vel_s[:, :, 3] * W  # rad/s
    omega_body_y = ball_vel_s[:, :, 4] * W
    omega_body_z = ball_vel_s[:, :, 5] * W
    spin_mag_rpm = np.sqrt(omega_body_x**2 + omega_body_y**2 + omega_body_z**2) * to_rpm
    true_roll_rpm = np.sqrt(omega_body_x**2 + omega_body_z**2) * to_rpm

    # ── Contact geometry ──
    sin_t = np.sin(theta_b)
    cos_t = np.cos(theta_b)

    ir_x_s = u_mat[:, ir_s:ir_s+1]       # shape (n,1)
    ir_y_s = u_mat[:, ir_s+1:ir_s+2]
    ir_z_s = u_mat[:, ir_s+2:ir_s+3]
    ir_tilt3 = u_mat[:, ir_s+3:ir_s+4]
    ir_tilt4 = u_mat[:, ir_s+4:ir_s+5]

    ir_radial_s = -ir_y_s * sin_t + ir_z_s * cos_t
    tilt_A1 = ir_tilt3 * cos_t + ir_tilt4 * sin_t

    # Inner groove center (nondim)
    r_gi0 = (d_m / 2 + (f_i - 0.5) * d_ball) / L
    x_gi0 = f_i * d_ball * np.sin(alpha_0) / L

    # Outer groove center
    r_go0 = (d_m / 2 + (f_o - 0.5) * d_ball) / L
    x_go0 = -f_o * d_ball * np.sin(alpha_0) / L  # negative for outer

    # Hertz stiffness (dimensional) — import from Julia or recalculate
    E_ball = 2.08e11; nu_ball = 0.3
    E_prime = E_ball / (1 - nu_ball**2)

    # Inner/outer curvature sums
    R_ball = d_ball / 2
    R_gi = f_i * d_ball
    R_go = f_o * d_ball
    R_inner_rolling = d_m / 2 - d_ball * np.cos(alpha_0) / 2
    R_outer_rolling = d_m / 2 + d_ball * np.cos(alpha_0) / 2

    sum_rho_i = 2/d_ball + 2/R_gi + 2/d_ball - 1/R_inner_rolling  # approximate
    sum_rho_o = 2/d_ball + 2/R_go + 2/d_ball - 1/R_outer_rolling

    # Use deformation-based Hertz: Q = Y * δ^1.5
    # Y = (4/3) * E' * √R_eff (for point contact)
    # For simplicity, let's use the nondim Hertz stiffness from kernel params
    # We'll approximate Y from the Julia printout
    Y_i = 3.2686e10   # N/m^1.5 (from Julia output)
    Y_o = 3.3999e10
    Y_i_star = Y_i * L**1.5 / Q_scale
    Y_o_star = Y_o * L**1.5 / Q_scale

    # drb distances
    drb_i = (f_i - 0.5) * d_ball / L
    drb_o = (f_o - 0.5) * d_ball / L

    # Contact geometry
    dx_i = x_gi0 + ir_x_s + r_gi0 / L * L * tilt_A1 / L * L / L  # This needs care...

    # Actually, simpler: just compute in nondim space matching kernel
    # dx_i = x_gi0★ + ir_x★ + r_gi0★ * tilt_A1 - x_b★
    # dr_i = r_gi0★ + dr_thermal★ + ir_radial★ - x_gi0★ * tilt_A1 - r_b★
    x_gi0_star = f_i * d_ball * np.sin(alpha_0) / L
    r_gi0_star = (d_m / 2 + (f_i - 0.5) * d_ball) / L

    dx_i = x_gi0_star + ir_x_s + r_gi0_star * tilt_A1 - x_b_s
    dr_i = r_gi0_star + ir_radial_s - x_gi0_star * tilt_A1 - r_b_s

    x_go0_star = -f_o * d_ball * np.sin(alpha_0) / L
    r_go0_star = (d_m / 2 + (f_o - 0.5) * d_ball) / L

    dx_o = x_go0_star - x_b_s
    dr_o = r_go0_star - r_b_s

    hypot_i = np.hypot(dx_i, dr_i)
    hypot_o = np.hypot(dx_o, dr_o)

    d_i_s = np.maximum(0.0, hypot_i - drb_i)
    d_o_s = np.maximum(0.0, hypot_o - drb_o)

    Q_inner = (Y_i_star * d_i_s**1.5) * Q_scale  # dimensional Newtons
    Q_outer = (Y_o_star * d_o_s**1.5) * Q_scale

    alpha_inner = np.degrees(np.arctan2(dx_i, dr_i))
    alpha_outer = np.degrees(np.arctan2(-dx_o, -dr_o))

    return {
        'time_s': t,
        'time_ms': t * 1e3,
        'ir_pos': ir_pos,
        'ir_vel': ir_vel,
        'cage_pos': cage_pos,
        'cage_vel': cage_vel,
        'cage_rpm': cage_rpm,
        'Q_inner': Q_inner,
        'Q_outer': Q_outer,
        'alpha_inner': alpha_inner,
        'alpha_outer': alpha_outer,
        'theta_deg': np.degrees(theta_b) % 360.0,
        'orbit_rpm': orbit_rpm,
        'spin_roll_rpm': true_roll_rpm,
        'spin_sideslip_rpm': omega_body_y * to_rpm,
        'spin_yaw_rpm': omega_body_x * to_rpm,
        'spin_pitch_rpm': omega_body_z * to_rpm,
        'spin_mag_rpm': spin_mag_rpm,
    }


def plot_race_displacement(data, out_dir):
    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    t_ms = data['time_ms']
    x_um = data['ir_pos'][:, 0] * 1e6
    y_um = data['ir_pos'][:, 1] * 1e6
    z_um = data['ir_pos'][:, 2] * 1e6

    axes[0].plot(t_ms, x_um, 'k-', lw=1.5)
    axes[0].set_title('Inner Race Axial Displacement (Vibration)', fontweight='bold')
    axes[0].set_ylabel('Axial X (μm)')

    axes[1].plot(t_ms, y_um, label='Radial Y', color='#1f77b4', lw=1.2, alpha=0.8)
    axes[1].plot(t_ms, z_um, label='Radial Z', color='#d62728', lw=1.2, alpha=0.8)
    axes[1].set_title('Inner Race Radial Displacement (Whirl & Vibration)', fontweight='bold')
    axes[1].set_xlabel('Time (ms)')
    axes[1].set_ylabel('Radial (μm)')
    axes[1].legend()

    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "01_Race_Displacement.png"), dpi=200)
    plt.close()
    print(f"  [OK] 01_Race_Displacement.png")


def plot_cage_dynamics(data, out_dir):
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    t_ms = data['time_ms']

    axes[0].plot(t_ms, data['cage_rpm'], label='Actual Cage Speed', color='orange', lw=2)
    axes[0].axhline(cage_theory_rpm, color='k', linestyle='--', label='Theoretical Epicyclic Speed')
    axes[0].set_title('Cage Speed vs Theoretical (Macro-Slip Monitor)', fontweight='bold')
    axes[0].set_xlabel('Time (ms)')
    axes[0].set_ylabel('Speed (RPM)')
    axes[0].legend()

    # Whirl orbit
    y_um = data['cage_pos'][:, 1] * 1e6
    z_um = data['cage_pos'][:, 2] * 1e6

    axes[1].plot(y_um, z_um, color='teal', alpha=0.7, lw=1.0)
    axes[1].plot(y_um[-1], z_um[-1], 'ro', label='End Point')

    # Approximate pilot clearance
    pilot_clr_um = 50  # approximate μm
    circle = plt.Circle((0, 0), pilot_clr_um, color='red', fill=False, linestyle='--', label='Pilot Clearance Limit')
    axes[1].add_artist(circle)
    limit = pilot_clr_um * 1.5
    axes[1].set_aspect('equal')
    axes[1].set_xlim(-limit, limit)
    axes[1].set_ylim(-limit, limit)
    axes[1].set_title('Cage Whirl Orbit (Y-Z Radial Plane)', fontweight='bold')
    axes[1].set_xlabel('Y Displacement (μm)')
    axes[1].set_ylabel('Z Displacement (μm)')
    axes[1].legend()

    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "02_Cage_Dynamics.png"), dpi=200)
    plt.close()
    print(f"  [OK] 02_Cage_Dynamics.png")


def plot_contact_loads(data, out_dir):
    fig = plt.figure(figsize=(14, 10))
    gs = GridSpec(2, 2, figure=fig)
    t_ms = data['time_ms']
    Q_i = data['Q_inner']
    Q_o = data['Q_outer']

    # Find heaviest loaded ball
    h_idx = np.argmax(np.mean(Q_o[-min(50, len(t_ms)):], axis=0))

    ax1 = fig.add_subplot(gs[0, :])
    ax1.plot(t_ms, Q_i[:, h_idx], label=f'Inner Race (Ball {h_idx})', lw=1.5, color='#1f77b4')
    ax1.plot(t_ms, Q_o[:, h_idx], label=f'Outer Race (Ball {h_idx})', lw=1.5, color='#d62728', linestyle='--')
    ax1.set_title('Transient Contact Loads (Tracking Heaviest Ball)', fontweight='bold')
    ax1.set_xlabel('Time (ms)')
    ax1.set_ylabel('Contact Load (N)')
    ax1.legend()

    ax2 = fig.add_subplot(gs[1, 0])
    balls = np.arange(Z)
    w = 0.35
    ax2.bar(balls - w/2, Q_i[-1, :], w, label='Inner', color='#1f77b4')
    ax2.bar(balls + w/2, Q_o[-1, :], w, label='Outer', color='#d62728')
    ax2.set_title('Load Distribution Across Balls (Final Step)', fontweight='bold')
    ax2.set_xlabel('Ball Number')
    ax2.set_ylabel('Load (N)')
    ax2.legend()

    ax3 = fig.add_subplot(gs[1, 1], polar=True)
    t_fin = np.radians(data['theta_deg'][-1, :])
    Q_fin = Q_o[-1, :]
    idx = np.argsort(t_fin)
    t_s = np.append(t_fin[idx], t_fin[idx[0]])
    Q_s = np.append(Q_fin[idx], Q_fin[idx[0]])
    ax3.plot(t_s, Q_s, color='darkred', marker='o', markersize=5)
    ax3.fill(t_s, Q_s, color='red', alpha=0.2)
    ax3.set_title('Outer Race Load Zone (Polar)', fontweight='bold', pad=15)

    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "03_Contact_Loads.png"), dpi=200)
    plt.close()
    print(f"  [OK] 03_Contact_Loads.png")


def plot_ball_kinematics(data, out_dir):
    fig, axes = plt.subplots(4, 1, figsize=(10, 16), sharex=True)
    t_ms = data['time_ms']

    mean_orbit = np.mean(data['orbit_rpm'], axis=1)
    axes[0].plot(t_ms, mean_orbit, 'k-', label='Mean Orbital Speed')
    axes[0].set_title('Ball Orbital Kinematics', fontweight='bold')
    axes[0].set_ylabel('Orbit Speed (RPM)')
    axes[0].legend()

    axes[1].plot(t_ms, data['spin_roll_rpm'][:, 0], color='#2ca02c', lw=1.2, label='Ball 0 True Roll (X-Z Vector)')
    axes[1].set_title('Ball True Rolling Speed (Healthy Kinematics)', fontweight='bold')
    axes[1].set_ylabel('Roll (RPM)')
    axes[1].legend()

    axes[2].plot(t_ms, data['spin_sideslip_rpm'][:, 0], color='#d62728', lw=1.2, label='Ball 0 Transverse Side-Slip')
    axes[2].set_title('Ball Gyroscopic Side-Slip & Drilling (Wear Indicator)', fontweight='bold')
    axes[2].set_ylabel('Speed (RPM)')
    axes[2].legend()

    axes[3].plot(t_ms, data['alpha_inner'][:, 0], label='Inner Contact Angle', color='#1f77b4')
    axes[3].plot(t_ms, data['alpha_outer'][:, 0], label='Outer Contact Angle', color='#d62728', linestyle='--')
    axes[3].axhline(alpha_0_deg, color='k', linestyle=':', label='Nominal Free Angle')
    axes[3].set_title('Dynamic Contact Angles (Ball 0) — Inner vs Outer', fontweight='bold')
    axes[3].set_xlabel('Time (ms)')
    axes[3].set_ylabel('Angle (Degrees)')
    axes[3].legend()

    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "04_Ball_Kinematics.png"), dpi=200)
    plt.close()
    print(f"  [OK] 04_Ball_Kinematics.png")


def print_bearing_frequencies():
    f_shaft = speed_rpm / 60.0
    f_c = f_shaft / 2 * (1 - gamma_cos)
    f_r = (1 / (2 * gamma)) * f_shaft * (1 - gamma_cos**2)
    f_ip = Z / 2 * f_shaft * (1 + gamma_cos)
    f_ep = Z / 2 * f_shaft * (1 - gamma_cos)

    print(f"=== Bearing Characteristic Frequencies ({speed_rpm:.0f} RPM) ===")
    print(f"  {'Frequency':<35} {'ADORE':>10} {'Unit':>6}")
    print(f"  {'-'*55}")
    print(f"  {'Inner ring, f_i':<35} {f_shaft:>10.3f} {'Hz':>6}")
    print(f"  {'Cage assembly, f_c':<35} {f_c:>10.3f} {'Hz':>6}")
    print(f"  {'Ball spin, f_r':<35} {f_r:>10.3f} {'Hz':>6}")
    print(f"  {'BPFI, f_ip':<35} {f_ip:>10.3f} {'Hz':>6}")
    print(f"  {'BPFO, f_ep':<35} {f_ep:>10.3f} {'Hz':>6}")
    print()


if __name__ == '__main__':
    npy_dir = r'e:\antigravity_folder\ADORE_V2_jl\julia_npy'

    if not os.path.exists(npy_dir):
        print(f"ERROR: {npy_dir} not found. Run Julia export first:")
        print(f"  julia --project=. export_30ms.jl")
        sys.exit(1)

    print("=" * 60)
    print("  Julia ADORE — Postprocessing Results")
    print("=" * 60)
    print(f"  Bearing: 7210B ({n_balls} balls)")
    print(f"  Load: Fa={F_axial:.0f} N, n={speed_rpm:.0f} RPM")
    print()

    print_bearing_frequencies()

    t, u_mat, L, V, W, Q = load_julia_npy(npy_dir)
    print(f"Loaded {len(t)} time steps from Julia .npy files")
    print(f"  Scales: L={L:.6f} m, V={V:.1f} m/s, W={W:.1f} rad/s, Q={Q:.1f} N")
    print()

    print("=== Extracting Kinematics ===")
    data = extract_kinematics(t, u_mat, float(L), float(V), float(W), float(Q))

    out_dir = r'e:\antigravity_folder\ADORE_V2_jl\results'
    os.makedirs(out_dir, exist_ok=True)

    print(f"\n=== Generating Plots → {out_dir} ===")
    plt.style.use('bmh')
    plot_race_displacement(data, out_dir)
    plot_cage_dynamics(data, out_dir)
    plot_contact_loads(data, out_dir)
    plot_ball_kinematics(data, out_dir)

    print(f"\n  Theory cage speed: {cage_theory_rpm:.1f} RPM")
    print(f"  Final cage speed:  {data['cage_rpm'][-1]:.1f} RPM ({data['cage_rpm'][-1]/cage_theory_rpm*100:.1f}%)")
    print("\nDone!")

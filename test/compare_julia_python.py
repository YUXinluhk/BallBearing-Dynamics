"""
compare_julia_python.py — Run Python 30ms sim with 7210B params (matching Julia),
load Julia CSV results, and generate comparison plots.
"""
import sys
sys.path.insert(0, r'e:\antigravity_folder')

import numpy as np
import math
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# ── ADORE Python imports ──
from ADORE_V2.bearing import BearingGeometry, MaterialParams
from ADORE_V2.ehl_traction import LubricantParams, traction_params_default
from ADORE_V2.cage import CageGeometry
from ADORE_V2.dynamics import BearingDynamics, SimulationConfig
from ADORE_V2.integrator import IntegratorConfig
from ADORE_V2.churning import ChurningParams

# ══════════════════════════════════════════════════════════════
# 1. Define 7210B bearing (MATCHING Julia config)
# ══════════════════════════════════════════════════════════════
def bearing_7210B():
    """7210B angular contact ball bearing — matching Julia params."""
    geom = BearingGeometry(
        d=12.7e-3,            # ball diameter [m]
        n_balls=16,           # 16 balls
        f_i=0.52,             # inner conformity
        f_o=0.53,             # outer conformity (Julia uses 0.53)
        d_m=70.0e-3,          # pitch diameter [m]
        alpha_0=np.radians(40.0),
        P_d=0.0,
        rho_ball=7800.0,
    )
    mat = MaterialParams(E=2.08e11, nu=0.3)
    return geom, mat


def lubricant_mil_l_23699():
    """MIL-L-23699 lubricant — matching Julia."""
    return LubricantParams(
        mu_0=0.01,
        alpha_pv=1.8e-8,
        beta_temp=0.03,
        T_0=373.0,
        K_th=0.14,
        rho_lub=860.0,
        c_p=2000.0,
    )


# ══════════════════════════════════════════════════════════════
# 2. Run Python simulation with matching parameters
# ══════════════════════════════════════════════════════════════
print("="*60)
print("  Python — 30ms Dynamic Simulation (7210B)")
print("="*60)

geom, mat = bearing_7210B()
lub = lubricant_mil_l_23699()
trac = traction_params_default()
cage_geom = CageGeometry.from_bearing(geom)

n_rpm = 10000.0
F_a = 2000.0
omega = n_rpm * 2 * np.pi / 60.0

print(f"  Bearing: 7210B ({geom.n_balls} balls)")
print(f"  Load: Fa={F_a:.0f} N, n={n_rpm:.0f} RPM")
print(f"  Duration: 30 ms")

dyn = BearingDynamics(geom, mat, lub, trac, cage_geom)

config = SimulationConfig(
    t_end=0.030,             # 30ms
    dt_output=50e-6,
    inner_race_speed=omega,
    F_axial=F_a,
    F_radial=0.0,
    delta_r_thermal=0,
    t_ramp_end=0.005,        # 5ms cosine ramp (Python default)
    zeta=0.10,
    integrator=IntegratorConfig(rtol=1e-4, atol=1e-7, h_max=1e-4),
    churning=ChurningParams(),
)

print("\nRunning Python simulation...")
result, scales = dyn.run(config, verbose=True)
print(f"\n  Python ODE success: {result.success}")
print(f"  Time points: {len(result.t)}")

# ── Extract Python time histories ──
t_py = result.t
n_t_py = len(t_py)
Z_py = geom.n_balls

# Scales
L_scale = scales.L
V_scale = scales.V
W_scale = scales.W

# State layout
n_pos_py = dyn.n_pos
ir_start = dyn._ir_start
ball_start_py = dyn._ball_start
cage_start_py = dyn._cage_start
n_ball_dof_py = dyn.n_ball_dof

# Inner race position (x, y, z): first 3 position DOFs
x_ir_py = result.y[ir_start, :] * L_scale        # axial [m]
y_ir_py = result.y[ir_start + 1, :] * L_scale
z_ir_py = result.y[ir_start + 2, :] * L_scale

# Inner race velocity
vx_ir_py = result.y[n_pos_py + ir_start, :] * V_scale
vy_ir_py = result.y[n_pos_py + ir_start + 1, :] * V_scale
vz_ir_py = result.y[n_pos_py + ir_start + 2, :] * V_scale

# Cage speed (θ̇_cage)
cage_omega_py = result.y[n_pos_py + cage_start_py + 3, :] * W_scale  # rad/s

# Ball orbital speed (θ̇_ball for each ball, numerical diff)
ball_orbital_py = np.zeros((n_t_py, Z_py))
for j in range(Z_py):
    theta_idx = ball_start_py + j * n_ball_dof_py + 2
    # θ is in nondim position, need velocity from velocity section
    theta_dot = result.y[n_pos_py + theta_idx, :] * W_scale
    ball_orbital_py[:, j] = theta_dot
mean_ball_orbital_py = np.mean(ball_orbital_py, axis=1)

# Ball 1 spin speed |ω|
ball1_spin_py = np.zeros(n_t_py)
for i in range(n_t_py):
    omega_idx_base = ball_start_py + 0 * n_ball_dof_py + 3  # ball 0, ωx position section
    wx = result.y[omega_idx_base, i] * W_scale
    wy = result.y[omega_idx_base + 1, i] * W_scale
    wz = result.y[omega_idx_base + 2, i] * W_scale
    ball1_spin_py[i] = np.sqrt(wx**2 + wy**2 + wz**2)

# Theory
gamma = geom.d * np.cos(geom.alpha_0) / geom.d_m
cage_theory_rads = 0.5 * omega * (1 - gamma)
cage_theory_rpm = cage_theory_rads * 30 / np.pi

print(f"\n  Python final cage: {cage_omega_py[-1]*30/np.pi:.1f} RPM (theory: {cage_theory_rpm:.1f})")
print(f"  Python final x_IR: {x_ir_py[-1]*1e6:.1f} μm")

# ══════════════════════════════════════════════════════════════
# 3. Load Julia CSV results
# ══════════════════════════════════════════════════════════════
julia_csv = r'e:\antigravity_folder\ADORE_V2_jl\julia_30ms_results.csv'
if not os.path.exists(julia_csv):
    print(f"\n  Julia CSV not found at {julia_csv}")
    print("  Run 'julia --project=. export_30ms.jl' first!")
    sys.exit(1)

print(f"\n  Loading Julia results from {julia_csv}")
jl_data = np.genfromtxt(julia_csv, delimiter=',', skip_header=1)
t_jl = jl_data[:, 0]
x_ir_jl = jl_data[:, 1]       # already in meters
y_ir_jl = jl_data[:, 2]
z_ir_jl = jl_data[:, 3]
vx_ir_jl = jl_data[:, 4]
cage_omega_jl = jl_data[:, 5]  # rad/s
mean_ball_orbital_jl = jl_data[:, 6]  # rad/s
ball1_spin_jl = jl_data[:, 7]  # rad/s

print(f"  Julia points: {len(t_jl)}")
print(f"  Julia final cage: {cage_omega_jl[-1]*30/np.pi:.1f} RPM")
print(f"  Julia final x_IR: {x_ir_jl[-1]*1e6:.1f} μm")

# ══════════════════════════════════════════════════════════════
# 4. Create comparison plots
# ══════════════════════════════════════════════════════════════
print("\n  Generating comparison figure...")

fig = plt.figure(figsize=(16, 14))
fig.suptitle(
    "ADORE — Julia vs Python Dynamic Simulation Comparison\n"
    f"7210B: {geom.n_balls} balls, {n_rpm:.0f} RPM, Fa={F_a:.0f} N, 30 ms",
    fontsize=14, fontweight='bold', y=0.98
)
gs = GridSpec(3, 2, hspace=0.35, wspace=0.30, top=0.91, bottom=0.06, left=0.08, right=0.95)

col_py = '#1f77b4'   # blue
col_jl = '#e74c3c'   # red
lw = 1.0

# ── Panel (a): Inner Race Axial Position ──
ax1 = fig.add_subplot(gs[0, 0])
ax1.plot(t_py * 1e3, x_ir_py * 1e6, color=col_py, lw=lw, alpha=0.8, label='Python')
ax1.plot(t_jl * 1e3, x_ir_jl * 1e6, color=col_jl, lw=lw, alpha=0.8, label='Julia', linestyle='--')
ax1.set_xlabel('Time [ms]')
ax1.set_ylabel('x_IR [μm]')
ax1.set_title('(a) Inner Race Axial Position')
ax1.legend(loc='best', fontsize=9)
ax1.grid(True, alpha=0.3)

# ── Panel (b): Inner Race Axial Velocity ──
ax2 = fig.add_subplot(gs[0, 1])
ax2.plot(t_py * 1e3, vx_ir_py, color=col_py, lw=lw, alpha=0.8, label='Python')
ax2.plot(t_jl * 1e3, vx_ir_jl, color=col_jl, lw=lw, alpha=0.8, label='Julia', linestyle='--')
ax2.set_xlabel('Time [ms]')
ax2.set_ylabel('v_IR_x [m/s]')
ax2.set_title('(b) Inner Race Axial Velocity')
ax2.legend(loc='best', fontsize=9)
ax2.grid(True, alpha=0.3)

# ── Panel (c): Cage Speed ──
ax3 = fig.add_subplot(gs[1, 0])
ax3.plot(t_py * 1e3, cage_omega_py * 30 / np.pi, color=col_py, lw=lw, alpha=0.8, label='Python')
ax3.plot(t_jl * 1e3, cage_omega_jl * 30 / np.pi, color=col_jl, lw=lw, alpha=0.8, label='Julia', linestyle='--')
ax3.axhline(cage_theory_rpm, color='green', linestyle=':', lw=1, alpha=0.7, label=f'Theory ({cage_theory_rpm:.0f} RPM)')
ax3.set_xlabel('Time [ms]')
ax3.set_ylabel('Cage Speed [RPM]')
ax3.set_title('(c) Cage Speed')
ax3.legend(loc='best', fontsize=9)
ax3.grid(True, alpha=0.3)

# ── Panel (d): Mean Ball Orbital Speed ──
ax4 = fig.add_subplot(gs[1, 1])
ax4.plot(t_py * 1e3, mean_ball_orbital_py * 30 / np.pi, color=col_py, lw=lw, alpha=0.8, label='Python')
ax4.plot(t_jl * 1e3, mean_ball_orbital_jl * 30 / np.pi, color=col_jl, lw=lw, alpha=0.8, label='Julia', linestyle='--')
ax4.axhline(cage_theory_rpm, color='green', linestyle=':', lw=1, alpha=0.7, label='Theory')
ax4.set_xlabel('Time [ms]')
ax4.set_ylabel('Mean Ball Orbital [RPM]')
ax4.set_title('(d) Mean Ball Orbital Speed')
ax4.legend(loc='best', fontsize=9)
ax4.grid(True, alpha=0.3)

# ── Panel (e): Ball 1 Spin Speed ──
ax5 = fig.add_subplot(gs[2, 0])
ax5.plot(t_py * 1e3, ball1_spin_py * 30 / np.pi, color=col_py, lw=lw, alpha=0.8, label='Python')
ax5.plot(t_jl * 1e3, ball1_spin_jl * 30 / np.pi, color=col_jl, lw=lw, alpha=0.8, label='Julia', linestyle='--')
ax5.set_xlabel('Time [ms]')
ax5.set_ylabel('Ball 1 |ω| [RPM]')
ax5.set_title('(e) Ball 1 Spin Speed')
ax5.legend(loc='best', fontsize=9)
ax5.grid(True, alpha=0.3)

# ── Panel (f): Summary text ──
ax6 = fig.add_subplot(gs[2, 1])
ax6.axis('off')
py_cage_final = cage_omega_py[-1] * 30 / np.pi
jl_cage_final = cage_omega_jl[-1] * 30 / np.pi

summary = (
    "Final State (t = 30 ms)\n"
    "─────────────────────────────────────\n"
    f"{'Metric':<28} {'Python':>10} {'Julia':>10}\n"
    f"{'─'*48}\n"
    f"{'x_IR [μm]':<28} {x_ir_py[-1]*1e6:>+10.1f} {x_ir_jl[-1]*1e6:>+10.1f}\n"
    f"{'v_IR [m/s]':<28} {vx_ir_py[-1]:>10.2e} {vx_ir_jl[-1]:>10.2e}\n"
    f"{'Cage [RPM]':<28} {py_cage_final:>10.1f} {jl_cage_final:>10.1f}\n"
    f"{'Cage/Theory':<28} {py_cage_final/cage_theory_rpm:>10.1%} {jl_cage_final/cage_theory_rpm:>10.1%}\n"
    f"{'Ball1 Spin [RPM]':<28} {ball1_spin_py[-1]*30/np.pi:>10.0f} {ball1_spin_jl[-1]*30/np.pi:>10.0f}\n"
    f"{'─'*48}\n"
    f"{'Theory cage':<28} {cage_theory_rpm:>10.0f} RPM\n"
)
ax6.text(0.05, 0.95, summary, fontfamily='monospace', fontsize=10,
         va='top', ha='left', transform=ax6.transAxes,
         bbox=dict(boxstyle='round,pad=0.5', facecolor='lightyellow', alpha=0.9))

out_path = r'e:\antigravity_folder\ADORE_V2_jl\julia_vs_python_comparison.png'
fig.savefig(out_path, dpi=150, bbox_inches='tight')
print(f"\n  Saved: {out_path}")
plt.close()

print("\n" + "="*60)
print("  COMPARISON COMPLETE")
print("="*60)

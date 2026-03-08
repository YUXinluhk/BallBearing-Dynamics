"""
export_python_30ms.py — Run Python 30ms sim with 7210B params (matching Julia) and export to CSV.
"""
import sys, os
sys.path.insert(0, r'e:\antigravity_folder')

import numpy as np
import math

from ADORE_V2.bearing import BearingGeometry, MaterialParams
from ADORE_V2.ehl_traction import LubricantParams, traction_params_default
from ADORE_V2.cage import CageGeometry
from ADORE_V2.dynamics import BearingDynamics, SimulationConfig
from ADORE_V2.integrator import IntegratorConfig
from ADORE_V2.churning import ChurningParams

# ── 7210B bearing (matching Julia) ──
geom = BearingGeometry(
    d=12.7e-3, n_balls=16, f_i=0.52, f_o=0.53,
    d_m=70.0e-3, alpha_0=np.radians(40.0), P_d=0.0, rho_ball=7800.0)
mat = MaterialParams(E=2.08e11, nu=0.3)
lub = LubricantParams(
    mu_0=0.01, alpha_pv=1.8e-8, beta_temp=0.03, T_0=373.0,
    K_th=0.14, rho_lub=860.0, c_p=2000.0)
trac = traction_params_default()
cage_geom = CageGeometry.from_bearing(geom)

n_rpm = 10000.0; F_a = 2000.0
omega = n_rpm * 2 * np.pi / 60.0

dyn = BearingDynamics(geom, mat, lub, trac, cage_geom)
config = SimulationConfig(
    t_end=0.030, dt_output=50e-6, inner_race_speed=omega,
    F_axial=F_a, F_radial=0.0, delta_r_thermal=0,
    t_ramp_end=0.005, zeta=0.10,
    integrator=IntegratorConfig(rtol=1e-4, atol=1e-7, h_max=1e-4),
    churning=ChurningParams())

print("Running Python 30ms (7210B, 16 balls, 10000 RPM, 2000N)...")
result, scales = dyn.run(config, verbose=True)
print(f"Success={result.success}, points={len(result.t)}")

# ── Extract data ──
# IMPORTANT: ADORE result.y has shape (n_t, n_state) — rows=time, cols=DOFs
t = result.t  # dimensional seconds
nP = dyn.n_pos
ir = dyn._ir_start
bs = dyn._ball_start
cs = dyn._cage_start
nd = dyn.n_ball_dof
Z = geom.n_balls
L, V, W = scales.L, scales.V, scales.W

x_ir = result.y[:, ir] * L
vx_ir = result.y[:, nP + ir] * V
cage_om = result.y[:, nP + cs + 3] * W  # cage theta_dot

# Ball orbital: use velocity DOF (theta_dot)
ball_orbit = np.zeros(len(t))
for j in range(Z):
    ball_orbit += result.y[:, nP + bs + j*nd + 2] * W
ball_orbit /= Z

# Ball 1 spin |omega|
b1s = np.zeros(len(t))
for i in range(len(t)):
    wx = result.y[i, bs + 3] * W  # omega in position section
    wy = result.y[i, bs + 4] * W
    wz = result.y[i, bs + 5] * W
    b1s[i] = np.sqrt(wx**2 + wy**2 + wz**2)

# Theory
gamma = geom.d * np.cos(geom.alpha_0) / geom.d_m
cage_th = 0.5 * omega * (1 - gamma)
print(f"Final cage: {cage_om[-1]*30/np.pi:.1f} RPM (theory: {cage_th*30/np.pi:.1f})")
print(f"Final x_IR: {x_ir[-1]*1e6:.1f} um")

# ── Export ──
out = np.column_stack([t, x_ir, np.zeros_like(t), np.zeros_like(t), vx_ir, cage_om, ball_orbit, b1s])
header = "t_s,x_ir_m,y_ir_m,z_ir_m,vx_ir_ms,cage_omega_rads,mean_ball_orbital_rads,ball1_spin_rads"
outf = r"e:\antigravity_folder\ADORE_V2_jl\python_30ms_results.csv"
np.savetxt(outf, out, delimiter=',', header=header, comments='')
print(f"Exported {len(t)} rows to {outf}")

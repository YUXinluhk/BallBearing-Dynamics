"""
plot_comparison.py — Create 6-panel Julia vs Python comparison figure from CSV data.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# ── Load data ──
jl = np.genfromtxt(r'e:\antigravity_folder\ADORE_V2_jl\julia_30ms_results.csv', delimiter=',', skip_header=1)
py = np.genfromtxt(r'e:\antigravity_folder\ADORE_V2_jl\python_30ms_results.csv', delimiter=',', skip_header=1)

t_jl, x_jl, vx_jl, cage_jl, orb_jl, spin_jl = jl[:,0], jl[:,1], jl[:,4], jl[:,5], jl[:,6], jl[:,7]
t_py, x_py, vx_py, cage_py, orb_py, spin_py = py[:,0], py[:,1], py[:,4], py[:,5], py[:,6], py[:,7]

# Theory (7210B: 16 balls, 10000 RPM, alpha=40 deg)
omega = 10000 * np.pi / 30
gamma = 12.7e-3 * np.cos(np.radians(40)) / 70e-3
cage_th = 0.5 * omega * (1 - gamma) * 30 / np.pi  # RPM

# ── Figure ──
fig = plt.figure(figsize=(16, 13))
fig.suptitle(
    'ADORE — Julia vs Python Dynamic Simulation\n'
    '7210B · 16 balls · 10 000 RPM · 2 000 N axial · 30 ms',
    fontsize=14, fontweight='bold', y=0.98)
gs = GridSpec(3, 2, hspace=0.38, wspace=0.28, top=0.91, bottom=0.06, left=0.08, right=0.95)

C_PY, C_JL = '#2563eb', '#dc2626'  # blue, red
LW = 1.2

# (a) Inner Race Axial Position
ax = fig.add_subplot(gs[0, 0])
ax.plot(t_py*1e3, x_py*1e6, C_PY, lw=LW, label='Python')
ax.plot(t_jl*1e3, x_jl*1e6, C_JL, lw=LW, label='Julia', ls='--')
ax.set(xlabel='Time [ms]', ylabel='x_IR [μm]', title='(a) Inner Race Axial Position')
ax.legend(fontsize=9); ax.grid(True, alpha=0.3)

# (b) Inner Race Axial Velocity
ax = fig.add_subplot(gs[0, 1])
ax.plot(t_py*1e3, vx_py, C_PY, lw=LW, label='Python')
ax.plot(t_jl*1e3, vx_jl, C_JL, lw=LW, label='Julia', ls='--')
ax.set(xlabel='Time [ms]', ylabel='v_x [m/s]', title='(b) Inner Race Axial Velocity')
ax.legend(fontsize=9); ax.grid(True, alpha=0.3)

# (c) Cage Speed
ax = fig.add_subplot(gs[1, 0])
ax.plot(t_py*1e3, cage_py*30/np.pi, C_PY, lw=LW, label='Python')
ax.plot(t_jl*1e3, cage_jl*30/np.pi, C_JL, lw=LW, label='Julia', ls='--')
ax.axhline(cage_th, color='#16a34a', ls=':', lw=1.5, label=f'Theory ({cage_th:.0f} RPM)')
ax.set(xlabel='Time [ms]', ylabel='Cage Speed [RPM]', title='(c) Cage Speed')
ax.legend(fontsize=9); ax.grid(True, alpha=0.3)

# (d) Mean Ball Orbital Speed
ax = fig.add_subplot(gs[1, 1])
ax.plot(t_py*1e3, orb_py*30/np.pi, C_PY, lw=LW, label='Python')
ax.plot(t_jl*1e3, orb_jl*30/np.pi, C_JL, lw=LW, label='Julia', ls='--')
ax.axhline(cage_th, color='#16a34a', ls=':', lw=1.5, label='Theory')
ax.set(xlabel='Time [ms]', ylabel='Mean Ball Orbital [RPM]', title='(d) Mean Ball Orbital Speed')
ax.legend(fontsize=9); ax.grid(True, alpha=0.3)

# (e) Ball 1 Spin Speed
ax = fig.add_subplot(gs[2, 0])
ax.plot(t_py*1e3, spin_py*30/np.pi, C_PY, lw=LW, label='Python')
ax.plot(t_jl*1e3, spin_jl*30/np.pi, C_JL, lw=LW, label='Julia', ls='--')
ax.set(xlabel='Time [ms]', ylabel='Ball 1 |ω| [RPM]', title='(e) Ball 1 Spin Speed')
ax.legend(fontsize=9); ax.grid(True, alpha=0.3)

# (f) Summary Table
ax = fig.add_subplot(gs[2, 1])
ax.axis('off')
py_c = cage_py[-1]*30/np.pi
jl_c = cage_jl[-1]*30/np.pi
txt = (
    "Final State (t = 30 ms)\n"
    "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
    f"{'Metric':<26} {'Python':>10} {'Julia':>10}\n"
    f"{'─'*46}\n"
    f"{'x_IR [μm]':<26} {x_py[-1]*1e6:>+10.1f} {x_jl[-1]*1e6:>+10.1f}\n"
    f"{'v_IR [m/s]':<26} {vx_py[-1]:>10.2e} {vx_jl[-1]:>10.2e}\n"
    f"{'Cage [RPM]':<26} {py_c:>10.1f} {jl_c:>10.1f}\n"
    f"{'Cage/Theory':<26} {py_c/cage_th:>10.1%} {jl_c/cage_th:>10.1%}\n"
    f"{'Ball1 |ω| [RPM]':<26} {spin_py[-1]*30/np.pi:>10.0f} {spin_jl[-1]*30/np.pi:>10.0f}\n"
    f"{'─'*46}\n"
    f"{'Theory cage':<26} {cage_th:>10.0f} RPM\n"
    f"\nNotes:\n"
    f"  Python: 5ms cosine ramp\n"
    f"  Julia:  no ramp (QS ICs at full speed)"
)
ax.text(0.05, 0.95, txt, fontfamily='monospace', fontsize=9.5,
        va='top', transform=ax.transAxes,
        bbox=dict(boxstyle='round,pad=0.5', fc='#fefce8', ec='#ca8a04', alpha=0.9))

out = r'e:\antigravity_folder\ADORE_V2_jl\julia_vs_python_comparison.png'
fig.savefig(out, dpi=150, bbox_inches='tight')
plt.close()
print(f'Saved: {out}')

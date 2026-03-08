"""
Force diagnostic: print forces on ball 1 at t=0 from Python kernel.
"""
import sys; sys.path.insert(0, r'e:\antigravity_folder')
import numpy as np
from ADORE_V2.bearing import BearingGeometry, MaterialParams
from ADORE_V2.ehl_traction import LubricantParams, traction_params_default
from ADORE_V2.cage import CageGeometry
from ADORE_V2.dynamics import BearingDynamics, SimulationConfig
from ADORE_V2.integrator import IntegratorConfig
from ADORE_V2.churning import ChurningParams

geom = BearingGeometry(d=12.7e-3,n_balls=16,f_i=0.52,f_o=0.53,d_m=70e-3,alpha_0=np.radians(40),P_d=0,rho_ball=7800)
mat = MaterialParams(E=2.08e11,nu=0.3)
lub = LubricantParams(mu_0=0.01,alpha_pv=1.8e-8,beta_temp=0.03,T_0=373,K_th=0.14,rho_lub=860,c_p=2000)
trac = traction_params_default()
cage = CageGeometry.from_bearing(geom)
dyn = BearingDynamics(geom,mat,lub,trac,cage)
c = SimulationConfig(t_end=0.001,dt_output=10e-6,inner_race_speed=10000*2*np.pi/60,F_axial=2000,F_radial=0,delta_r_thermal=0,t_ramp_end=0.005,zeta=0.10,integrator=IntegratorConfig(rtol=1e-4,atol=1e-7,h_max=1e-4),churning=ChurningParams())
r,s = dyn.run(c,verbose=False)

# Print state evolution at very early points
nP = dyn.n_pos; cs = dyn._cage_start; bs = dyn._ball_start; nd = dyn.n_ball_dof
omega = 10000*2*np.pi/60
gamma = geom.d*np.cos(geom.alpha_0)/geom.d_m
cage_th = 0.5*omega*(1-gamma)

print("Python early time evolution:")
print(f"{'t_ms':>8} {'cage_RPM':>10} {'cage%':>7} {'ball_td_RPM':>12} {'ball_vx':>10} {'x_ir_um':>8}")
for i in range(min(20, len(r.t))):
    t = r.t[i]
    cage_w = r.y[i, nP+cs+3] * s.W
    b_td = r.y[i, nP+bs+2] * s.W
    b_vx = r.y[i, nP+bs+0] * s.V
    x_ir = r.y[i, 0] * s.L
    print(f"{t*1e3:8.4f} {cage_w*30/np.pi:10.1f} {cage_w/cage_th*100:7.1f}% {b_td*30/np.pi:12.1f} {b_vx:10.2e} {x_ir*1e6:8.2f}")

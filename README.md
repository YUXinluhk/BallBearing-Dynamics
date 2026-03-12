# ADORE.jl — Six-Degree-of-Freedom Transient Thermo-Elasto-Hydrodynamic Bearing Dynamics Solver

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

**ADORE** (Advanced Dynamic Orbiting Rolling Element) is a high-fidelity transient bearing dynamics solver implemented in Julia. It integrates six-degree-of-freedom rigid body mechanics, tribological contact micro-mechanics, and thermal network modelling into a single implicit ODE system with closed-loop energy conservation.

## Key Features

| Module | Description |
|--------|-------------|
| **6-DOF Dynamics** | Quaternion-based 3D rotation for each ball with gyroscopic precession, coupled to translational DOFs and cage dynamics |
| **84-Point TEHD Contact** | Sub-contact resolution of Heathcote creep, spin micro-slip, Bair–Winer limiting shear stress, and Archard flash temperature on an 11×11 grid |
| **4-Node LPTN Thermal Network** | Lumped-parameter thermal network (inner race, outer race, balls, lubricant) with dynamic thermal expansion |
| **Energy Audit** | ODE-integrated energy conservation verified to < 0.02% error |
| **Sparse Jacobian** | Graph-colored finite differencing with analytically derived sparsity pattern (83% sparse); AD-ready kernel |

## Validated Against

- **Palmgren** friction torque model
- **SKF** catalog speed limits  
- **Harris** analytical cage speed theory (< 0.2% error)
- **NASA TP-2275** experimental thermal data (35mm & 120mm bearings, DN up to 3×10⁶)
- Multi-bearing validation: 7008C, 7010C, 7014C, 7210B

## Quick Start

```julia
# 1. Clone and activate
using Pkg
Pkg.activate(".")
Pkg.instantiate()

# 2. Run a single case
include("run_single.jl")

# 3. Run transient plots (spin nutation, cage whirl, thermal settling)
include("run_transient_plots.jl")
```

## Project Structure

```
├── src/
│   ├── Config/        # TOML parser and type definitions
│   ├── Dynamics/      # ODE kernel, driver, state vector, Jacobian sparsity
│   ├── Contact/       # Hertzian geometry, EHL film, traction integrals
│   └── ADORE.jl       # Module entry point
├── inputs/            # TOML configuration files for various bearing cases
├── paper/             # Manuscript (LaTeX)
├── run_single.jl      # Standard simulation + 14 diagnostic plots
└── run_transient_plots.jl  # Extended transient simulation (2s) + 3 publication figures
```

## Citation

If you use ADORE in your research, please cite:

> Lu, Y.X. *ADORE: A High-Fidelity Six-Degree-of-Freedom Transient Thermo-Elasto-Hydrodynamic Bearing Dynamics Solver with Closed-Loop Energy Conservation.* (2025).

## License

MIT

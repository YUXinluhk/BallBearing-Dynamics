# Peer Review for Applied Sciences

## Manuscript
`ADORE: A High-Fidelity Six-Degree-of-Freedom Transient Thermo-Elasto-Hydrodynamic Bearing Dynamics Solver with Closed-Loop Energy Conservation`

## Abstract
This manuscript presents a high-fidelity transient dynamic simulation framework for angular contact ball bearings. The main contributions include a six-degree-of-freedom rigid-body model for each rolling element with quaternion kinematics, a thermo-elasto-hydrodynamic contact traction model, an 84-point sub-contact integration scheme, thermal-network coupling, and implicit stiff ODE integration using FBDF. The principal conclusions are that the solver reproduces cage speed and energy conservation well for a 7210B bearing and shows good agreement with NASA TP-2275 heat-generation data for the 35 mm bearing case.

Overall, the topic is important and the technical ambition is high. The manuscript has the potential to become a valuable methodology/software paper in bearing dynamics. However, in its current form, the validation claims are stronger than the evidence presented, and several key modeling and reproducibility issues remain insufficiently resolved.

## Review Decision
**Major Revision**

## Major Comments

### 1. Validation claims are stronger than the evidence currently supports
The manuscript claims strong agreement for both 35 mm and 120 mm NASA bearing cases, but the cross-bearing prediction for the 120 mm case shows a systematic underprediction of about 70 percent. This is not a minor discrepancy; it indicates that the present model does not yet demonstrate robust cross-scale predictive capability.

The abstract and conclusions should therefore be revised to distinguish clearly between:

- successful agreement on the calibrated 35 mm case, and
- unsuccessful blind extrapolation to the 120 mm case.

If the authors wish to retain a generalization claim, they should provide additional independent validation beyond the currently presented dataset.

### 2. Calibration and validation are not sufficiently separated
The 35 mm NASA case appears to be used, at least in part, for parameter calibration of traction-related quantities. If so, the resulting low error should not be presented as an independent validation in the strict sense. In addition, the NASA data are extracted from published figures, which introduces additional uncertainty.

The manuscript should explicitly separate:

- parameter calibration,
- model verification/self-consistency checks, and
- independent validation.

At minimum, the authors should hold out some operating points as a validation set not used for parameter determination. They should also state the uncertainty introduced by digitization of NASA curves.

### 3. Internal inconsistency in the thermal model description reduces technical credibility
The paper alternates between describing the thermal model as a four-node and a five-node network. This inconsistency appears in the introduction, body text, figure captions, and conclusions. Because the thermal network size directly affects the state dimension, energy balance interpretation, and reproducibility, this issue must be corrected carefully throughout the manuscript.

Similarly, the cage state definition should be checked for consistency, since angular velocity appears to be mixed into generalized coordinates in one place.

The authors should perform a thorough consistency audit of:

- thermal network size,
- state vector definition,
- state-count formula,
- figure captions,
- abstract, and
- conclusions.

### 4. The treatment of internal clearance and preload is not sufficiently defined
For angular contact ball bearings, internal clearance or preload is a central modeling ingredient. The manuscript notes that thermal expansion can compress internal clearance and may trigger thermal runaway, but the initial clearance/preload state and its role in contact geometry are not defined clearly enough.

This is a major mechanical modeling issue. Without this information, it is difficult to judge whether the load distribution, contact angle evolution, and thermo-mechanical coupling are physically reliable.

The paper should explicitly report:

- initial radial and/or axial clearance,
- any assembly preload,
- boundary conditions for inner and outer rings,
- how clearance/preload enters the contact geometry equations, and
- a sensitivity analysis showing how predicted loads, heat generation, and temperature respond to clearance/preload variation.

### 5. The reproducibility claim is not yet fully supported by the parameter disclosure
The manuscript states that the framework is described in sufficient detail for independent reproduction, but several quantities in the governing equations are not fully linked to the tabulated parameter set. For example, the traction formulation introduces parameters such as `A`, `B`, `C`, `D`, and `Lambda_LSS`, while the appendix tables list `kappa_0`, `kappa_infty`, `kappa_m`, and `mu_spin` without a clear mapping.

To make the work genuinely reproducible, the manuscript should provide:

- a complete mapping from all symbolic parameters in the equations to the numerical values used in the implementation,
- units and sources for each parameter,
- justification for damping, conductance, and regularization terms,
- ODE solver tolerances and convergence criteria in full detail, and
- sensitivity of the results to numerical settings such as tolerances and maximum time step.

### 6. The temperature comparison is weakened by the use of thermally accelerated capacitances
The manuscript reports temperature comparison after reducing the thermal capacitances by a factor of 100 in order to reach a quasi-steady condition within a 150 ms simulation window. This procedure changes the thermal time scale and therefore cannot be treated as a direct one-to-one validation of measured temperature evolution.

This part should be reframed more cautiously, for example as a pseudo-steady extrapolation rather than strict thermal validation. If the authors want to preserve a strong temperature-validation claim, they should either:

- simulate the full physical thermal transient using real capacitances, or
- restrict the main validation claim to heat generation and treat temperature ordering as a qualitative secondary result.

### 7. Several presented checks are useful, but they are primarily self-consistency checks rather than external validation
Comparison against Harris cage speed, energy conservation, and Palmgren/SKF friction estimates is helpful and should remain in the paper. However, these are mainly verification or consistency checks rather than strong external validation of bearing dynamics.

For a manuscript positioned as a high-fidelity bearing-dynamics solver, the paper would be substantially strengthened by adding at least one of the following:

- direct comparison with experimentally observed cage motion or whirl behavior,
- comparison with measured temperature histories rather than endpoint values only,
- comparison with friction torque measurements, or
- quantitative comparison of dynamic signatures such as orbit behavior, load redistribution, or spin-axis evolution from published experiments.

## Minor Comments

### Technical Writing
- The manuscript still uses the `Tribology International` journal field and should be fully adapted if this version is intended for `Applied Sciences`.
- The author affiliation and email currently appear to be placeholders and should be finalized.
- Terminology should be standardized throughout the paper, especially regarding `four-node` versus `five-node` thermal network and related state definitions.
- Some claims are written too categorically and should be softened where the supporting evidence is partial.

### Methodology Presentation
- The numerical workflow should be easier to reproduce by adding a compact algorithmic summary of initialization, transient integration, Jacobian construction, and post-processing.
- The origin of the NASA data points extracted from figures should be described more transparently, including any digitization procedure and uncertainty.
- Zero-speed cases should be explained more carefully in the energy table rather than simply marked with dashes.

### Figures and Tables
- The current figures read more like internal diagnostic plots than publication-oriented validation figures.
- The manuscript would benefit from more direct experiment-versus-simulation overlays with clear legends, normalized axes where appropriate, and uncertainty indication.
- The table definitions for energy-balance quantities should be clarified so the reader immediately understands whether values refer to absolute mismatch or relative error.

### Structure
- The logical flow is generally reasonable, but parts of the Results and Discussion sections overlap.
- The manuscript would read more clearly if the authors separated:
  - what was modeled,
  - what was numerically verified,
  - what was experimentally validated, and
  - what physical insights can actually be claimed from the current evidence.

## Conclusion
This is a promising manuscript with substantial technical depth and a relevant contribution to the field of bearing dynamics. The integration of rigid-body rolling-element dynamics, TEHD traction, thermal coupling, and energy auditing into a unified solver is a meaningful step forward.

However, the present manuscript requires major revision before it can reach the standard expected for acceptance in `Applied Sciences`. The most important actions are to reduce overstatement in the validation claims, clearly separate calibration from validation, fully specify clearance/preload treatment, correct internal inconsistencies in the thermal model description, and strengthen the external validation narrative. If these issues are addressed carefully, the work could become a strong and useful contribution.

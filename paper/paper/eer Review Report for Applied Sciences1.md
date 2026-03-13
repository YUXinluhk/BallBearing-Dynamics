# Peer Review Report for Applied Sciences

## Manuscript Information

**Title**: ADORE: A High-Fidelity Six-Degree-of-Freedom Transient Thermo-Elasto-Hydrodynamic Bearing Dynamics Solver with Closed-Loop Energy Conservation

**Field**: Bearing Dynamics / Tribology

**Journal Target**: Applied Sciences

---

## Abstract

This manuscript presents ADORE, an open-source high-fidelity transient simulation framework for angular contact ball bearings. The work combines a six-degree-of-freedom dynamic model for each rolling element, quaternion-based rotational tracking, a lumped thermal network, an 84-point thermo-elasto-hydrodynamic contact integral, and an implicit stiff ODE solver with sparse Jacobian acceleration. The main conclusions are that the framework reproduces cage kinematics with high accuracy relative to Harris theory, maintains excellent internal energy-balance consistency, and predicts NASA thermal data for the 35-mm bearing within a small mean error.

Overall, the manuscript is ambitious, technically sophisticated, and potentially impactful. It demonstrates substantial modeling depth and a commendable effort toward reproducibility. However, several important issues remain regarding internal consistency, validation scope, parameter identifiability, and journal-format suitability. These issues should be addressed before the manuscript can be considered for acceptance.

---

## Review Decision

**Major Revision**

---

## Major Comments

### 1. Internal inconsistency in the thermal-network description

The manuscript contains a clear inconsistency regarding the number of thermal nodes.

- The abstract describes a **four-node** lumped-parameter thermal network.
- The thermal-model section also defines four temperature states: inner race, outer race, ball, and oil.
- However, the Introduction and Conclusions refer to a **five-node** thermal network.

This inconsistency must be corrected throughout the manuscript. If the ambient is only a boundary condition and not a dynamic state, the model is four-node and should be described consistently as such. If a fifth node exists in the implementation, then the mathematical formulation must be revised accordingly.

### 2. Cross-bearing prediction for the 120-mm bearing reveals a major limitation

The cross-bearing validation shows a systematic underprediction of heat generation by roughly 70% for the 120-mm NASA bearing. This is not a small discrepancy; it indicates that the present model does not yet generalize reliably across bearing size using fixed traction and churning parameters.

The discussion currently offers plausible qualitative explanations, but the analysis remains too speculative. The manuscript would be significantly strengthened by at least one quantitative sensitivity test showing how much of the discrepancy can be attributed to:

- traction-parameter variation,
- churning/fill-fraction assumptions,
- thermal boundary assumptions.

At minimum, this section should be reframed more explicitly as a **model limitation** rather than a successful predictive extension.

### 3. The “first-principles” claim should be moderated

The manuscript repeatedly emphasizes first-principles prediction, yet key traction parameters appear to be selected differently for the 7210B and NASA 35-mm cases. This creates tension between the modeling claim and the actual calibration practice.

Please clarify:

- how the traction parameters were chosen,
- whether they were fitted to experiments,
- whether they come from rheological measurements or literature,
- and which parts of the framework are genuinely predictive versus parameter-tuned.

A more precise formulation would be to state that the solver provides **first-principles sub-contact resolution within a semi-empirical constitutive framework**, unless a stronger parameter-identification basis is supplied.

### 4. Validation remains limited in independence

The current validation set combines:

- Harris theory for cage speed,
- internal energy-balance closure,
- NASA TP-2275 thermal data,
- Palmgren and SKF empirical comparisons.

This is useful, but still limited. Harris theory is not experimental validation, and energy conservation is an internal consistency check rather than external evidence of physical correctness. The manuscript would benefit from at least one additional independent validation source, such as:

- published experimental temperature/torque data from another study,
- comparison with a known bearing-dynamics code,
- or the already cited `Takabi and Khonsari (2013)` dataset if appropriate.

### 5. Lack of uncertainty or sensitivity analysis

For a numerical paper of this complexity, the absence of sensitivity analysis is a significant weakness. The results likely depend strongly on several difficult-to-identify parameters, including:

- thermal conductances,
- traction coefficients,
- oil fill fraction,
- smoothing parameter for contact regularization,
- quaternion stabilization gain.

The manuscript should include at least a compact sensitivity study showing whether the main reported conclusions remain stable under reasonable parameter variation.

### 6. Journal formatting is not aligned with Applied Sciences

The manuscript is presently formatted using an `elsarticle` structure and explicitly names `Tribology International` as the journal. This is incompatible with submission to Applied Sciences. The manuscript must be reformatted using the proper Applied Sciences/MDPI template before final submission.

---

## Minor Comments

### 1. Thermal-network terminology should be made fully consistent

The manuscript alternates between “four-node” and “five-node” descriptions. This appears in more than one section and should be corrected globally.

### 2. State-vector presentation could be clarified

The state-dimension count is mathematically consistent, but the mapping between the listed generalized coordinates and the final dimension formula is not immediately transparent. Consider presenting a short table summarizing the DOF and state contributions from:

- inner race,
- each ball,
- cage,
- thermal states,
- energy-audit states.

This would improve readability.

### 3. Clarify the relationship between traction coefficient and shear stress

The manuscript introduces an isothermal traction coefficient and later uses shear stress expressions, but the connection between them is not written explicitly enough. Please state the constitutive link more clearly so that the TEHD formulation can be reproduced without ambiguity.

### 4. Figures need stronger self-contained captions

Several figures appear technically relevant but are not fully self-explanatory from the captions alone. The captions should indicate more clearly:

- what quantity is plotted,
- whether the curves are transient or steady-state,
- and what physical conclusion the reader should draw.

This is especially important for the NASA comparison figures and the Case 4 diagnostics.

### 5. Qualitative comparison with Gupta should be labeled carefully

The section comparing the results with Gupta’s ADORE is interesting but remains mostly qualitative. Since no direct code-to-code or case-to-case quantitative benchmark is available, the wording should avoid implying stronger validation than is actually demonstrated.

### 6. Provide rationale for the Baumgarte stabilization parameter

The choice of the quaternion stabilization gain is explained qualitatively, but the manuscript would benefit from brief evidence that the selected value does not distort the dynamics or introduce artificial numerical damping.

### 7. Add missing or underused references

`Takabi and Khonsari (2013)` appears in the bibliography but is not discussed in the manuscript. Either use this reference substantively or remove it. Also, the manuscript would benefit from a few more recent references from the last five years to situate the work in the current literature.

### 8. Churning-model citation should be completed

The discussion mentions the `Goksem-Hargreaves` churning model, but the corresponding reference is missing. Please add the full citation.

### 9. Code availability statement should be strengthened

The code-availability section is appreciated, but it would be better to include:

- a version identifier or commit hash,
- a permanent archive DOI if available,
- and a brief note on how the manuscript cases can be reproduced from the repository.

### 10. Numerical-method choice could be motivated more explicitly

The manuscript uses FBDF, which is a reasonable choice for a stiff multiscale system. Still, the paper would benefit from a short explanation of why this method was preferred over alternative stiff integrators such as implicit Runge-Kutta or Rosenbrock-type methods.

---

## Comments on Scientific Rigor

The manuscript generally shows strong scientific maturity, especially in its attention to multiscale stiffness, coordinate representation, and thermal-mechanical coupling. The treatment of quaternion kinematics, contact smoothing, and vectorial pocket friction is thoughtful and technically strong.

However, the central scientific concern is not whether the model is sophisticated, but whether its predictive capability is robust outside the cases where parameters were effectively adapted. The manuscript would be more convincing if the authors distinguish more carefully between:

- validated physics,
- calibrated constitutive assumptions,
- and unresolved scale effects.

---

## Comments on Methodology

The methodology is described in unusually good detail compared with many simulation papers, and the appendix helps reproducibility. The explanation of the Jacobian sparsity strategy and stiff integration approach is a clear strength.

That said, reproducibility would be improved further by:

- explicitly documenting how traction parameters are selected,
- clarifying all thermal conductance assumptions,
- and stating the exact workflow used to generate each validation figure and table.

---

## Comments on Data Interpretation

The interpretations of cage-speed agreement, energy-balance closure, and thermal ordering are generally reasonable and grounded in bearing theory. The discussion of co-rotating entrainment velocity and Heathcote conformity is particularly strong.

The main interpretive weakness is that the manuscript occasionally overstates what the data prove. For example:

- agreement with Harris theory confirms kinematic consistency, but not full dynamic validity;
- energy conservation confirms internal numerical consistency, but not physical correctness;
- successful fit to one NASA bearing size does not establish cross-scale generality.

These distinctions should be made more explicit.

---

## Comments on Technical Writing

The paper is generally well written, technically precise, and more polished than most first submissions in this field. The author clearly understands the underlying mechanics and tribology.

Still, several improvements would help:

- use fully consistent terminology for the thermal network,
- moderate a few claims that sound stronger than the evidence supports,
- and ensure all sections align with the target journal style and formatting.

---

## Comments on Figures, Tables, and Structure

The overall structure of the paper is logical and effective: governing equations, numerical implementation, validation, and interpretation follow a sensible progression. The tables are useful and mostly well designed.

The figures would benefit from more informative captions and stronger integration into the narrative. In particular, the manuscript should help the reader understand why each figure is included and what exact physical feature is being validated or illustrated.

---

## Conclusion

This is a substantial and technically impressive manuscript with clear potential for publication after careful revision. The coupled dynamic-thermal-lubrication framework is ambitious and the implementation appears to be of high quality. The work is especially strong in model formulation, numerical treatment, and reproducibility intent.

Before the manuscript can reach the standard expected for acceptance in Applied Sciences, the authors should address the internal inconsistencies, clarify the role of calibrated parameters, strengthen the validation framing, and more honestly delimit the current predictive scope of the model. With these revisions, the paper could become a valuable contribution to the bearing-dynamics literature.

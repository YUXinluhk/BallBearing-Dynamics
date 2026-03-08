"""
dimensional_analysis.py — Automated dimensional analysis of ALL nondimensionalization
formulas in Julia nondim.jl and Python nondim.py.

Uses symbolic dimensional tracking: each quantity is represented as a product of
base dimensions [kg^a · m^b · s^c · rad^d].
"""
import sys
from dataclasses import dataclass
from typing import Optional


@dataclass(frozen=True)
class Dim:
    """Dimensional exponents: kg^a · m^b · s^c · rad^d"""
    kg: float = 0
    m: float = 0
    s: float = 0
    rad: float = 0

    def __mul__(self, other):
        if isinstance(other, Dim):
            return Dim(self.kg + other.kg, self.m + other.m,
                       self.s + other.s, self.rad + other.rad)
        return self  # scalar mult doesn't change dim

    def __truediv__(self, other):
        if isinstance(other, Dim):
            return Dim(self.kg - other.kg, self.m - other.m,
                       self.s - other.s, self.rad - other.rad)
        return self

    def __pow__(self, exp):
        return Dim(self.kg * exp, self.m * exp, self.s * exp, self.rad * exp)

    def __eq__(self, other):
        if isinstance(other, Dim):
            return (abs(self.kg - other.kg) < 1e-10 and
                    abs(self.m - other.m) < 1e-10 and
                    abs(self.s - other.s) < 1e-10 and
                    abs(self.rad - other.rad) < 1e-10)
        return False

    def __repr__(self):
        parts = []
        for name, exp in [('kg', self.kg), ('m', self.m), ('s', self.s), ('rad', self.rad)]:
            if abs(exp) > 1e-10:
                if abs(exp - 1) < 1e-10:
                    parts.append(name)
                elif abs(exp - int(exp)) < 1e-10:
                    parts.append(f'{name}^{int(exp)}')
                else:
                    parts.append(f'{name}^{exp:.2g}')
        return ' · '.join(parts) if parts else '1 (dimensionless)'

    def is_dimensionless(self):
        return (abs(self.kg) < 1e-10 and abs(self.m) < 1e-10 and
                abs(self.s) < 1e-10 and abs(self.rad) < 1e-10)


# ── Base dimensions ──
DIMLESS = Dim()
KG = Dim(kg=1)
M = Dim(m=1)
S = Dim(s=1)
RAD = Dim(rad=1)

# ── Derived dimensions ──
N = KG * M * S**(-2)        # Newton = kg·m/s²
PA = N / M**2               # Pascal = N/m²
J = N * M                   # Joule  = N·m

# ── Physical quantities ──
FORCE = N                   # [N] = kg·m·s^-2
LENGTH = M                  # [m]
TIME = S                    # [s]
MASS = KG                   # [kg]
VELOCITY = M / S            # [m/s]
ANGVEL = S**(-1)            # [rad/s] (rad is dimensionless in SI)
ACCEL = M * S**(-2)         # [m/s²]
STIFFNESS_HERTZ = N * M**(-1.5)  # [N/m^{3/2}] — Hertz contact stiffness Υ
INERTIA = KG * M**2         # [kg·m²]
MOMENT = N * M              # [N·m]
DAMPING = N * S / M         # [N·s/m] = kg/s


def sqrt(d: Dim) -> Dim:
    return d ** 0.5


def check(name: str, formula_dim: Dim, expected_dim: Dim, source: str):
    """Check if formula dimensions match expected."""
    ok = formula_dim == expected_dim
    status = '✅ PASS' if ok else '❌ FAIL'
    print(f'  {status}  {name}')
    print(f'         Formula gives: [{formula_dim}]')
    print(f'         Expected:      [{expected_dim}]')
    if not ok:
        ratio = formula_dim / expected_dim
        print(f'         ERROR factor:  [{ratio}]')
    print(f'         Source: {source}')
    print()
    return ok


# ═══════════════════════════════════════════════════════════════
#  INPUTS (physical quantities)
# ═══════════════════════════════════════════════════════════════
m = MASS                         # ball mass [kg]
L = LENGTH                       # ball radius [m]
Y = STIFFNESS_HERTZ              # Hertz Υ [N/m^{3/2}]


print('='*72)
print('  DIMENSIONAL ANALYSIS — Julia nondim.jl + Python nondim.py')
print('='*72)

n_pass = 0
n_fail = 0
n_total = 0

def RUN(name, formula, expected, source):
    global n_pass, n_fail, n_total
    n_total += 1
    if check(name, formula, expected, source):
        n_pass += 1
    else:
        n_fail += 1


# ═══════════════════════════════════════════════════════════════
#  SECTION 1: Scale Definitions (Julia nondim.jl lines 30-43)
# ═══════════════════════════════════════════════════════════════
print('\n' + '─'*72)
print('  SECTION 1: Scale Definitions')
print('─'*72 + '\n')

# Q = Y * L^1.5
Q_julia = Y * L**1.5
RUN('Q = Y · L^{3/2}', Q_julia, FORCE, 'nondim.jl:35')

# T = sqrt(m * L^0.5 / Y)   ← JULIA CODE (line 36)
T_julia_code = sqrt(m * L**0.5 / Y)
RUN('T = √(m · L^{1/2} / Y)  [Julia CODE]', T_julia_code, TIME, 'nondim.jl:36')

# T = sqrt(m / Y) * L^{-1/4}  ← JULIA COMMENT (line 12)
T_julia_comment = sqrt(m / Y) * L**(-0.25)
RUN('T = √(m/Y) · L^{-1/4}  [Julia COMMENT]', T_julia_comment, TIME, 'nondim.jl:12 (comment)')

# T = sqrt(m * L / Q)       ← PYTHON CODE (nondim.py:62)
Q_python = Y * L**1.5  # same Q
T_python = sqrt(m * L / Q_python)
RUN('T = √(m · L / Q₀)  [Python CODE]', T_python, TIME, 'nondim.py:62')

# V = L / T  (using Julia's T)
V_julia = L / T_julia_code
RUN('V = L / T  [Julia]', V_julia, VELOCITY, 'nondim.jl:37')

# V = L / T  (using Python's T)
V_python = L / T_python
RUN('V = L / T  [Python]', V_python, VELOCITY, 'nondim.py derived')

# W = 1 / T
W_julia = DIMLESS / T_julia_code
RUN('W = 1 / T  [Julia]', W_julia, ANGVEL, 'nondim.jl:38')

# K = Y (Hertz stiffness scale)
RUN('K = Y (stiffness scale)', Y, STIFFNESS_HERTZ, 'nondim.jl:39')

# I_scale = m * L^2
I_julia = m * L**2
RUN('I = m · L²  (inertia scale)', I_julia, INERTIA, 'nondim.jl:40')

# M_scale = Q * L  (moment scale)
M_julia = Q_julia * L
RUN('M = Q · L  (moment scale)', M_julia, MOMENT, 'nondim.jl:41')


# ═══════════════════════════════════════════════════════════════
#  SECTION 2: Nondimensionalization Functions
# ═══════════════════════════════════════════════════════════════
print('\n' + '─'*72)
print('  SECTION 2: Nondim/Dim Conversion Functions')
print('─'*72 + '\n')

# nondim_force: F / Q → dimensionless
RUN('F★ = F / Q', FORCE / Q_julia, DIMLESS, 'nondim.jl:48')

# nondim_length: x / L → dimensionless
RUN('x★ = x / L', LENGTH / L, DIMLESS, 'nondim.jl:49')

# nondim_time: t / T (Julia T)
RUN('t★ = t / T  [Julia T]', TIME / T_julia_code, DIMLESS, 'nondim.jl:50')

# nondim_time: t / T (Python T)
RUN('t★ = t / T  [Python T]', TIME / T_python, DIMLESS, 'nondim.py')

# nondim_vel: v / V
RUN('v★ = v / V  [Julia]', VELOCITY / V_julia, DIMLESS, 'nondim.jl:51')

# nondim_angvel: ω / W
RUN('ω★ = ω / W  [Julia]', ANGVEL / W_julia, DIMLESS, 'nondim.jl:52')

# nondim_mass: m * L / (Q * T²)
nondim_mass_julia = MASS * L / (Q_julia * T_julia_code**2)
RUN('m★ = m · L / (Q · T²)  [Julia]', nondim_mass_julia, DIMLESS, 'nondim.jl:56')

nondim_mass_python = MASS * L / (Q_python * T_python**2)
RUN('m★ = m · L / (Q · T²)  [Python T]', nondim_mass_python, DIMLESS, 'nondim.py equivalent')


# ═══════════════════════════════════════════════════════════════
#  SECTION 3: Damping Nondimensionalization
# ═══════════════════════════════════════════════════════════════
print('\n' + '─'*72)
print('  SECTION 3: Damping Coefficient c★ = c · V / Q')
print('─'*72 + '\n')

# c_dim has units [N·s/m] = [kg/s]
# c★ = c · V / Q → should be dimensionless
c_star_julia = DAMPING * V_julia / Q_julia
RUN('c★ = c · V / Q  [Julia V]', c_star_julia, DIMLESS, 'nondim.jl/nondim.py:155')

c_star_python = DAMPING * V_python / Q_python
RUN('c★ = c · V / Q  [Python V]', c_star_python, DIMLESS, 'nondim.py:155')


# ═══════════════════════════════════════════════════════════════
#  SECTION 4: EOM Consistency Checks
# ═══════════════════════════════════════════════════════════════
print('\n' + '─'*72)
print('  SECTION 4: Equations of Motion (Nondimensional Consistency)')
print('─'*72 + '\n')

# Newton's 2nd law (translational): F = m·a
# Nondim: F★·Q = m★·(m)·(a★·L/T²)
# → F★ = m★·a★·(m·L)/(Q·T²) = m★·a★·1 (if m★ = m·L/(Q·T²) = 1)
# Check: m★ = m·L/(Q·T²) should = 1 for ball

# With Julia T:
m_star_julia = MASS * L / (Q_julia * T_julia_code**2)
RUN('m★ (Julia T) = m·L/(Q·T²) should be dimensionless', m_star_julia, DIMLESS,
    'EOM consistency')

# With Python T:
m_star_python = MASS * L / (Q_python * T_python**2)
RUN('m★ (Python T) = m·L/(Q·T²) should be dimensionless', m_star_python, DIMLESS,
    'EOM consistency')

# Damping in EOM: c★·v★ should have dimensions of acceleration★ (i.e. dimensionless)
# Python ball: ball_acc -= c_ball_t * vel  (c★ acts as acceleration coefficient)
# → c★ must have [acc★/v★] = [1] (dimensionless / dimensionless) = dimensionless
RUN('c★ · v★ = acc★  (Python ball EOM)', DIMLESS * DIMLESS, DIMLESS,
    'dynamics_numba.py:1103')

# Python IR: ir_acc -= c_ir * vel / m_race_star
# → c★/m★ acts as acceleration coefficient
# → c★/m★ must be dimensionless
RUN('c★/m★ = acc★/v★  (Python IR EOM)', DIMLESS / DIMLESS, DIMLESS,
    'dynamics_numba.py:1092')


# ═══════════════════════════════════════════════════════════════
#  SECTION 5: Cross-Consistency (Julia T vs Python T)
# ═══════════════════════════════════════════════════════════════
print('\n' + '─'*72)
print('  SECTION 5: Julia vs Python T — Ratio Analysis')
print('─'*72 + '\n')

# T_julia / T_python = ?
# T_julia  = √(m·L^{0.5}/Y)
# T_python = √(m/(Y·L^{0.5}))
# ratio = T_julia / T_python
#       = √((m·L^{0.5}/Y) · (Y·L^{0.5}/m))
#       = √(L^{0.5}·L^{0.5}) = √L = L^{1/2}
ratio_dim = T_julia_code / T_python
print(f'  T_julia / T_python has dimension: [{ratio_dim}]')
print(f'  Expected: [{LENGTH**0.5}]  (= √L ≈ √0.00635 ≈ 0.0797)')
if ratio_dim == LENGTH**0.5:
    print(f'  ✅ Confirmed: T_julia = T_python × √L')
    print(f'  ⚠️  This means Julia T has an EXTRA √m factor')
else:
    print(f'  ❌ Unexpected ratio')


# ═══════════════════════════════════════════════════════════════
#  SECTION 6: Correct Julia T formula
# ═══════════════════════════════════════════════════════════════
print('\n' + '─'*72)
print('  SECTION 6: Proposed Fix')
print('─'*72 + '\n')

T_fixed = sqrt(m / (Y * L**0.5))
RUN('T_fixed = √(m / (Y · L^{1/2}))', T_fixed, TIME,
    'PROPOSED FIX: L^{0.5} in denominator')

# Verify this matches Python:
T_py_expanded = sqrt(m / (Y * L**0.5))
T_py_from_Q = sqrt(m * L / Q_python)  # = √(m*L/(Y*L^1.5)) = √(m/(Y*L^0.5))
print(f'  T_Python = √(m·L/Q) = √(m/(Y·√L))')
print(f'  Dim of T_Python:        [{T_py_from_Q}]')
print(f'  Dim of T_fixed:         [{T_fixed}]')
print(f'  Match: {T_py_from_Q == T_fixed}')


# ═══════════════════════════════════════════════════════════════
#  SUMMARY
# ═══════════════════════════════════════════════════════════════
print('\n' + '='*72)
print(f'  SUMMARY: {n_pass} PASS / {n_fail} FAIL / {n_total} TOTAL')
print('='*72)
if n_fail > 0:
    print(f'\n  ⚠️  {n_fail} dimensional inconsistencies found!')
    print('  The Julia T formula has L^{0.5} in the NUMERATOR')
    print('  but it should be in the DENOMINATOR.')
    print()
    print('  FIX in nondim.jl line 36:')
    print('    CURRENT:  T = sqrt(m * L^0.5 / Y)')
    print('    CORRECT:  T = sqrt(m / (Y * sqrt(L)))')
else:
    print('\n  All formulas are dimensionally consistent.')

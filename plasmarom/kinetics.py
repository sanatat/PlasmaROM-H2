"""Dimensionless, illustrative H2-equivalent excitation/oxidation network.

Not a calibrated plasma chemistry, hydrogen combustion, or Boltzmann solver.
Species are conserved *hydrogen fuel equivalents*, not literal mole fractions.
  F0 -> F1 -> F2 -> R -> P (excited pools, reactive pool, converted pool)
Plasma drives F0 -> F1 and F1 -> F2. Relaxation and oxidation compete.
"""
from dataclasses import dataclass
import numpy as np
from scipy.integrate import solve_ivp, trapezoid

FULL_NAMES = ('H2(v0) proxy', 'H2(v1) proxy', 'H2(v2) proxy', 'Reactive pool', 'Converted pool')
RED_NAMES = ('H2(v0) proxy', 'H2(v1+v2) lumped pool', 'Reactive pool', 'Converted pool')

@dataclass(frozen=True)
class Parameters:
    excitation: float = 2.5
    second_excitation: float = 2.0
    relaxation10: float = 1.2
    relaxation21: float = 12.0
    activated_oxidation: float = 2.0
    background_oxidation: float = 0.03
    radical_oxidation: float = 1.3
    radical_quench: float = 0.07
    pulse_strength: float = 1.0
    pulse_width: float = 0.34


def pulse(t, p=Parameters()):
    """Smooth dimensionless electric excitation history; NOT a discharge model."""
    return p.pulse_strength * (np.exp(-((t - 1.2) / p.pulse_width)**2) +
                               0.75 * np.exp(-((t - 3.0) / p.pulse_width)**2))


def rates(t, p=Parameters(), spatial=1.0):
    drive = pulse(t, p) * spatial
    return p.excitation * drive, p.second_excitation * drive


def rhs_full(t, y, p=Parameters(), spatial=1.0):
    f0, f1, f2, radical, product = y
    a01, a12 = rates(t, p, spatial)
    q01 = a01 * f0
    q10 = p.relaxation10 * f1
    q12 = a12 * f1
    q21 = p.relaxation21 * f2
    q2r = p.activated_oxidation * f2
    q0r = p.background_oxidation * f0
    qrp = p.radical_oxidation * radical
    qr0 = p.radical_quench * radical
    return np.stack((-q01 + q10 - q0r + qr0,
                     q01 - q10 - q12 + q21,
                     q12 - q21 - q2r,
                     q2r + q0r - qrp - qr0,
                     qrp), axis=0)


def qssa_partition(t, p=Parameters(), spatial=1.0):
    """Algebraic fast-manifold fraction F2/(F1+F2), on a toy transfer network.

    Using a12*F1 = (relaxation21 + activated_oxidation)*F2,
    with B=F1+F2 gives F1=B*(1-frac), F2=B*frac.
    Not valid if the fast rate denominator is zero or timescale separation fails.
    """
    denominator = p.relaxation21 + p.activated_oxidation
    if denominator <= 0 or not np.isfinite(denominator):
        raise ValueError('QSSA requires a positive finite fast-state outflow rate')
    _, a12 = rates(t,p,spatial)
    if np.any(np.asarray(a12)<0):
        raise ValueError('QSSA excitation rate cannot be negative')
    return a12 / (denominator + a12)


def rhs_reduced(t, y, p=Parameters(), spatial=1.0):
    """Conservative fast-pool lumping: y=[F0,B=F1+F2,R,P].

    The original model's ``F1`` variable was incorrectly used both as F1
    and B. All F1/F2-specific transitions now consistently use the QSSA
    partition of B. The approximation is still poor when the drive changes
    faster than the F2 relaxation timescale; no global accuracy is claimed.
    """
    f0, b, radical, product = y
    a01, _ = rates(t, p, spatial)
    fraction2 = qssa_partition(t, p, spatial)
    f1 = (1.-fraction2)*b
    f2 = fraction2*b
    q01 = a01*f0
    q10 = p.relaxation10*f1
    q2r = p.activated_oxidation*f2
    q0r = p.background_oxidation*f0
    qrp = p.radical_oxidation*radical
    qr0 = p.radical_quench*radical
    return np.stack((-q01+q10-q0r+qr0,
                     q01-q10-q2r,
                     q2r+q0r-qrp-qr0,
                     qrp), axis=0)


def simulate(p=Parameters(), reduced=False, method='BDF', time_end=6.0, n_points=501, rtol=1e-8, atol=1e-11):
    y0 = [1., 0., 0., 0.] if reduced else [1., 0., 0., 0., 0.]
    fun = rhs_reduced if reduced else rhs_full
    return solve_ivp(lambda t,y: fun(t,y,p), (0, time_end), y0,
                     t_eval=np.linspace(0, time_end, n_points),
                     method=method, rtol=rtol, atol=atol)


def full_observables(solution):
    """Return [F0, B=F1+F2, R, P], matching the reduced state definitions."""
    return np.stack((solution.y[0],solution.y[1]+solution.y[2],
                     solution.y[3],solution.y[4]),axis=0)


def pca_variance(matrix):
    """Variance retained by k principal components of row=samples data."""
    centered = matrix - matrix.mean(axis=0, keepdims=True)
    _, s, vt = np.linalg.svd(centered, full_matrices=False)
    var = s*s
    ratio = var / var.sum() if var.sum() else np.zeros_like(var)
    return ratio, vt


def reaction_flux_screen(t, y, p=Parameters()):
    """Positive integrated reaction-pathway fluxes, NOT a DRGEP implementation.

    Screening ranks toy channels by integrated throughput. Back-and-forth
    excitation can rank highly despite little net progress toward P.
    """
    drive = pulse(t, p)
    f0,f1,f2,r,_ = y
    fluxes = {
      'F0 -> F1 (plasma)': p.excitation*drive*f0,
      'F1 -> F0 (relaxation)': p.relaxation10*f1,
      'F1 -> F2 (plasma)': p.second_excitation*drive*f1,
      'F2 -> F1 (relaxation)': p.relaxation21*f2,
      'F2 -> R (activated)': p.activated_oxidation*f2,
      'F0 -> R (background)': p.background_oxidation*f0,
      'R -> P (oxidation)': p.radical_oxidation*r,
      'R -> F0 (quench)': p.radical_quench*r
    }
    return sorted(((name,float(trapezoid(rate,t))) for name,rate in fluxes.items()),
                  key=lambda pair:pair[1],reverse=True)


def qssa_reconstruct(red, p=Parameters()):
    """Approximate F2 recovery from QSSA output, not independent simulation.

    Reduced F1 is regarded as a lumped (F1+F2) pool; split it according to
    algebraic QSSA F2/F1 = k12/(k21+k2r). Use ONLY for visualization and
    checking discrepancies; do not claim this proves physical validity.
    """
    if red.y.shape[0]!=4:
        raise ValueError('Expected 4-state reduced solution')
    fraction2=qssa_partition(red.t,p)
    f2=red.y[1]*fraction2
    f1=red.y[1]-f2
    return np.vstack((red.y[0],f1,f2,red.y[2],red.y[3]))

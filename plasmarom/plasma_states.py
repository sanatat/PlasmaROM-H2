"""Dimensioned three-level H2 vibrational kinetic example (v=0,1,2).

NIST H2 molecular spectroscopic constants are used for vibrational G(v).
Electron excitation uses tabulated cross sections integrated against an
ASSUMED Maxwellian electron EEDF. Not an EEDF solver or a discharge model.
Synthetic cross sections for tests/demo are fictional and explicitly labelled.
"""
from dataclasses import dataclass
from pathlib import Path
import csv
import json
import numpy as np
from scipy.integrate import solve_ivp, trapezoid
from scipy.constants import e as E_CHARGE, m_e, h, c

# NIST CCCBDB H2 vibrational constants: term value in inverse centimetres.
W_E_CM=4401.213
W_EX_E_CM=121.336
W_EY_E_CM=0.8129


def levels_eV():
    v=np.arange(3,dtype=float)
    x=v+0.5
    term=W_E_CM*x-W_EX_E_CM*x**2+W_EY_E_CM*x**3
    return (term-term[0])*h*c*100/E_CHARGE


def read_cross_section_csv(path):
    """Read energy_eV,sigma_m2 table; strictly ascending energy grid."""
    with open(path,newline='',encoding='utf-8') as f:
        rows=list(csv.DictReader(f))
    try:
        table=np.asarray([[float(row['energy_eV']),float(row['sigma_m2'])] for row in rows])
    except (KeyError,ValueError,TypeError) as exc:
        raise ValueError('CSV must contain numeric energy_eV and sigma_m2 columns') from exc
    if table.ndim!=2 or table.shape[0]<3 or table.shape[1]!=2:
        raise ValueError('Expected >=3 cross-section rows')
    if not np.all(np.isfinite(table)) or np.any(table<0) or not np.all(np.diff(table[:,0])>0) or table[0,0]!=0:
        raise ValueError('Require finite nonnegative, strictly ascending energies starting at zero')
    return table


def maxwell_rate(table, te_eV, threshold_eV, *, tail_tolerance=1e-4):
    """Return <sigma v> [m3/s] using a Maxwellian energy distribution.

    Rate = sqrt(8e/(pi*m_e))*Te^-3/2 * int sigma(E)*E*exp(-E/Te)dE;
    Te and E are in eV. Extrapolation ABOVE table max is not performed;
    user must cover the Maxwell tail for desired precision.
    """
    table=np.asarray(table,float)
    if (table.ndim!=2 or table.shape[1]!=2 or table.shape[0]<3 or
        not np.all(np.isfinite(table)) or np.any(table<0) or
        table[0,0]!=0 or not np.all(np.diff(table[:,0])>0) or
        not np.isfinite(te_eV) or te_eV<=0 or not np.isfinite(threshold_eV) or threshold_eV<0):
        raise ValueError('Invalid cross section, threshold or electron temperature')
    if not (0 < tail_tolerance < 1):
        raise ValueError('tail_tolerance must be in (0,1)')
    grid=table[:,0];sigma=table[:,1].copy()
    if threshold_eV >= grid[-1]:
        raise ValueError('Cross-section grid does not cover the excitation threshold')
    # The normalized Maxwellian energy-weight above Emax. This is a
    # distribution-tail coverage test, NOT an upper bound on rate truncation
    # when the unknown cross section beyond Emax could increase.
    x=grid[-1]/te_eV
    tail_mass=(1+x)*np.exp(-x)
    if tail_mass > tail_tolerance:
        raise ValueError(f'Insufficient EEDF energy coverage: normalized tail {tail_mass:.3g} > {tail_tolerance:g}; extend cross-section table')
    sigma[grid<threshold_eV]=0
    # Do not integrate nonzero interpolated values between the last below-
    # threshold grid node and the threshold itself.
    energies=np.unique(np.r_[grid,threshold_eV])
    cross=np.interp(energies,grid,sigma,left=0,right=0)
    cross[energies<threshold_eV]=0
    if threshold_eV>0:
        cross[energies==threshold_eV]=0
    integral=trapezoid(cross*energies*np.exp(-energies/te_eV),energies)
    return float(np.sqrt(8*E_CHARGE/(np.pi*m_e))*integral/te_eV**1.5)



def maxwell_reversible_rates(table, te_eV, threshold_eV, *, degeneracy_lower=1.,
                             degeneracy_upper=1., tail_tolerance=1e-4):
    """Pair of excitation/de-excitation rate coefficients for a Maxwellian bath.

    Both rates use ONE threshold-aware quadrature grid. The stable reverse
    integral is evaluated with exp(-(E-Delta)/Te); the forward rate is its
    Boltzmann-scaled counterpart. This enforces microscopic reversibility:
    k_down / k_up = (g_lower/g_upper) * exp(Delta/Te).

    Cross-section and Maxwellian hypotheses are externally imposed. This is
    not a Boltzmann EEDF solution, nor measured electron-impact data.
    """
    # This routine reuses maxwell_rate's validation, EEDF tail coverage guard
    # and input conventions; the returned forward quadrature is not reused.
    maxwell_rate(table, te_eV, threshold_eV, tail_tolerance=tail_tolerance)
    if (not np.isfinite(degeneracy_lower) or not np.isfinite(degeneracy_upper)
            or degeneracy_lower <= 0 or degeneracy_upper <= 0):
        raise ValueError('Both level degeneracies must be finite and positive')
    table=np.asarray(table, dtype=float)
    E0=table[:,0]
    threshold=float(threshold_eV)
    E=np.unique(np.r_[E0,threshold])
    sigma=np.interp(E,E0,table[:,1])
    sigma[E <= threshold]=0.0
    prefactor=np.sqrt(8*E_CHARGE/(np.pi*m_e))/float(te_eV)**1.5
    reverse_base=prefactor*trapezoid(sigma*E*np.exp(-(E-threshold)/te_eV),E)
    reverse=float(reverse_base*degeneracy_lower/degeneracy_upper)
    forward=float(reverse_base*np.exp(-threshold/te_eV))
    return forward,reverse

def synthetic_table(threshold,scale=1e-20):
    """Invented illustration ONLY. NEVER describe as H2 experimental data."""
    E=np.linspace(0,35,3501)
    sigma=np.zeros_like(E)
    active=E>threshold
    sigma[active]=scale*(1-np.exp(-(E[active]-threshold)/1.5))*np.exp(-(E[active]-threshold)/35)
    return np.column_stack((E,sigma))


@dataclass(frozen=True)
class PlasmaConfig:
    te_eV:float=2.
    ne_m3:float=1e17
    h2_m3:float=1e22
    quench10_s:float=6e3
    quench21_s:float=6e3
    duration_s:float=2e-4


def integrate_levels(table01,table12,config=PlasmaConfig()):
    pars=(config.te_eV,config.ne_m3,config.h2_m3,config.quench10_s,config.quench21_s,config.duration_s)
    if any(not np.isfinite(z) for z in pars) or min(pars[:3])<=0 or min(pars[3:5])<0 or pars[5]<=0:
        raise ValueError('Invalid plasma parameters')
    levels=levels_eV();de01=levels[1]-levels[0];de12=levels[2]-levels[1]
    k01,k10=maxwell_reversible_rates(table01,config.te_eV,de01)
    k12,k21=maxwell_reversible_rates(table12,config.te_eV,de12)
    up01=config.ne_m3*k01;up12=config.ne_m3*k12
    down10=config.ne_m3*k10;down21=config.ne_m3*k21
    q10=config.quench10_s;q21=config.quench21_s
    def rhs(t,y):
        n0,n1,n2=y[:3]
        f01=up01*n0;f10=down10*n1
        f12=up12*n1;f21=down21*n2
        q1=q10*n1;q2=q21*n2
        # cumulative energy per m3: electron reservoir change and gas heat.
        pe=-config.h2_m3*E_CHARGE*(de01*(f01-f10)+de12*(f12-f21))
        pg=config.h2_m3*E_CHARGE*(de01*q1+de12*q2)
        return [-f01+f10+q1, f01-f10-f12+f21-q1+q2,
                f12-f21-q2,pe,pg]
    times=np.linspace(0,config.duration_s,301)
    sol=solve_ivp(rhs,(0,config.duration_s),[1.,0.,0.,0.,0.],t_eval=times,
                  method='BDF',rtol=1e-10,atol=1e-13)
    if not sol.success:
        raise RuntimeError(sol.message)
    vibration=config.h2_m3*E_CHARGE*(levels[:,None]*sol.y[:3]).sum(axis=0)
    energy_residual=np.max(np.abs(vibration+sol.y[3]+sol.y[4]))
    energy_scale=max(1.,float(np.max(np.abs(vibration))),float(np.max(abs(sol.y[3]))))
    return {'times':sol.t,'fractions':sol.y[:3], 'electron_energy_Jm3':sol.y[3],
            'gas_heat_Jm3':sol.y[4],'levels_eV':levels,'excitation_rates_m3s':(k01,k12),
            'max_population_sum_error':float(np.max(abs(sol.y[:3].sum(axis=0)-1))),
            'min_population':float(sol.y[:3].min()),
            'relative_energy_closure':float(energy_residual/energy_scale)}


def run_plasma(out, path01=None,path12=None, source01=None,source12=None,
               config=PlasmaConfig()):
    synthetic=(path01 is None and path12 is None)
    if not synthetic and (not path01 or not path12 or not source01 or not source12):
        raise ValueError('Two CSVs plus both source/provenance references are mandatory')
    gaps=np.diff(levels_eV())
    a=synthetic_table(gaps[0]) if synthetic else read_cross_section_csv(path01)
    b=synthetic_table(gaps[1]) if synthetic else read_cross_section_csv(path12)
    results=integrate_levels(a,b,config)
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    label='SYNTHETIC CROSS SECTIONS; NOT PHYSICAL VALIDATION' if synthetic else 'USER-SUPPLIED CROSS SECTIONS; VALIDITY NOT VERIFIED'
    record={'data_status':label,'provenance':('invented demonstration' if synthetic else {'v0_v1':source01,'v1_v2':source12}),
            'levels_eV':results['levels_eV'].tolist(),
            'rate_coefficients_m3s':list(results['excitation_rates_m3s']),
            'max_population_sum_error':results['max_population_sum_error'],
            'min_population':results['min_population'],
            'relative_energy_closure':results['relative_energy_closure'],
            'limitations':'Assumed Maxwell EEDF, no E/N or Poisson equation, no electron balance or calibrated plasma, no ignition chemistry'}
    (out/'plasma_results.json').write_text(json.dumps(record,indent=2))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(7,4.5))
    for i in range(3):ax.plot(results['times']*1e3,results['fractions'][i],label=f'H2(v={i})')
    ax.legend();ax.set(xlabel='Time [ms]',ylabel='Vibrational population fraction',title=label)
    fig.tight_layout();fig.savefig(out/'plasma_states.png',dpi=160);plt.close(fig)
    return record

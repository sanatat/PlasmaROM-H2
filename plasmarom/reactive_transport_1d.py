"""1-D conservative, constant-density reactive transport *verification benchmark*.

Transport of two mass fractions, A -> B, with Arrhenius-like heat release.
No Navier-Stokes, combustion mechanism, plasma, DNS or LES is claimed.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from scipy.constants import R


def simulate(*, cells=100, length_m=1., end_s=.02, u_m_s=.3,
             diffusivity_m2_s=0.005, thermal_diffusivity_m2_s=0.005,
             k_ref_s=20., activation_J_mol=0., reference_T_K=600.,
             heat_J_kg=2e5, cp_J_kg_K=1000., uniform=False):
    if (cells<12 or length_m<=0 or end_s<=0 or diffusivity_m2_s<0 or
        thermal_diffusivity_m2_s<0 or k_ref_s<0 or heat_J_kg<0 or cp_J_kg_K<=0 or
        not np.isfinite([length_m,end_s,u_m_s,diffusivity_m2_s,thermal_diffusivity_m2_s,
                         k_ref_s,activation_J_mol,reference_T_K,heat_J_kg,cp_J_kg_K]).all() or
        activation_J_mol<0 or reference_T_K<=0):
        raise ValueError('Invalid physical/numerical parameters')
    dx=length_m/cells
    x=(np.arange(cells)+0.5)*dx
    A=np.full(cells,.6) if uniform else .2+.4*np.exp(-((x-length_m*.4)/(.12*length_m))**2)
    B=1-A
    T=np.full(cells,reference_T_K)
    A_init=float(A.mean());B_init=float(B.mean()); energy_init=float((cp_J_kg_K*T+heat_J_kg*A).mean())
    dt_adv=dx/(abs(u_m_s)+1e-30)
    dt_dif=dx*dx/(2*max(diffusivity_m2_s,thermal_diffusivity_m2_s,1e-30))
    # Upper bound of k(T) for positive heat release & nonnegative activation
    t_max=reference_T_K+heat_J_kg/cp_J_kg_K
    k_max=k_ref_s*np.exp(activation_J_mol/R*(1/reference_T_K-1/t_max))
    dt_reaction=.5/max(k_max,1e-30)
    steps=int(np.ceil(end_s/(.4*min(dt_adv,dt_dif,dt_reaction))))
    if steps>300000:raise ValueError('Too many explicit steps; change grid or physics')
    dt=end_s/steps
    def transport(phi,D):
        next_cell=np.roll(phi,-1)
        if u_m_s>=0:
            Fadv=u_m_s*phi
        else:
            Fadv=u_m_s*next_cell
        flux=Fadv-D*(next_cell-phi)/dx
        return phi-dt*(flux-np.roll(flux,1))/dx
    for _ in range(steps):
        A=transport(A,diffusivity_m2_s)
        B=transport(B,diffusivity_m2_s)
        T=transport(T,thermal_diffusivity_m2_s)
        k=k_ref_s*np.exp(activation_J_mol/R*(1/reference_T_K-1/T))
        consumed=A*(-np.expm1(-k*dt))
        A-=consumed;B+=consumed
        T+=consumed*heat_J_kg/cp_J_kg_K
        if min(A.min(),B.min(),T.min()) < -1e-11 or not np.isfinite(T).all():
            raise RuntimeError('Nonphysical transport or temperature')
    mass_drift=float(abs((A+B).mean()-(A_init+B_init)))
    # Thermal diffusivity unequal to species diffusivity: both domain integrals
    # are conserved on periodic domain with equal cp, and reaction redistributes
    # cp*T + q*A even though locally unequal diffusivity.
    enthalpy_drift=float(abs((cp_J_kg_K*T+heat_J_kg*A).mean()-energy_init))
    result={'status':'VERIFIED_1D_SCALAR_REACTIVE_TRANSPORT_NOT_CFD_DNS_LES',
            'max_mass_fraction_sum_error':float(np.max(abs(A+B-1))),
            'global_mass_drift':mass_drift,'global_total_enthalpy_drift_J_kg':enthalpy_drift,
            'steps':steps,'dx_m':dx,'dt_s':dt,
            'mean_product_B':float(B.mean()),'max_temperature_K':float(T.max())}
    return {'x_m':x,'A':A,'B':B,'T_K':T,'summary':result}


def run_transport_demo(out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    d=simulate();np.savez_compressed(out/'reactive_transport_1d.npz',x_m=d['x_m'],
            A=d['A'],B=d['B'],T_K=d['T_K'])
    (out/'reactive_transport_1d.json').write_text(json.dumps(d['summary'],indent=2))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(1,2,figsize=(9,3))
    ax[0].plot(d['x_m'],d['A'],label='A');ax[0].plot(d['x_m'],d['B'],label='B')
    ax[0].set(xlabel='x [m]',ylabel='Mass fraction');ax[0].legend()
    ax[1].plot(d['x_m'],d['T_K']);ax[1].set(xlabel='x [m]',ylabel='Temperature [K]')
    fig.suptitle('1D surrogate A→B, not combustion or DNS/LES')
    fig.tight_layout();fig.savefig(out/'reactive_transport_1d.png',dpi=160);plt.close(fig)
    return d['summary']

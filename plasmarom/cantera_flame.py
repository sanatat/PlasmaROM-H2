"""Optional CHEMICAL 1D premixed H2/air flame using Cantera FreeFlame.

This is a reference integration exercise, not an internally developed 1D
flow solver and not plasma-assisted combustion. Requires Cantera >= 3.2.
"""
import json
from pathlib import Path
import numpy as np
from .cantera_demo import _ct


def run_flame(out, mechanism='h2o2.yaml', width_m=0.03, initial_T_K=300.,
              pressure_atm=1.0, phi=1.0, refine_slope=0.08, refine_curve=0.14):
    if (not all(np.isfinite(x) and x>0 for x in (width_m,initial_T_K,pressure_atm,phi,refine_slope,refine_curve))
        or not mechanism):
        raise ValueError('Positive finite flame input parameters required')
    ct=_ct()
    gas=ct.Solution(mechanism)
    gas.set_equivalence_ratio(phi,'H2','O2:1,N2:3.76')
    gas.TP=float(initial_T_K),float(pressure_atm)*ct.one_atm
    f=ct.FreeFlame(gas,width=width_m)
    f.transport_model='mixture-averaged'
    f.set_refine_criteria(ratio=3,slope=refine_slope,curve=refine_curve)
    f.solve(loglevel=0,auto=True)
    speed=float(f.velocity[0]); Tmax=float(np.max(f.T))
    if not np.isfinite(speed) or speed<=0 or not np.isfinite(Tmax) or Tmax<=initial_T_K:
        raise RuntimeError('Nonphysical or unconverged flame result')
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    h2_idx=gas.species_index('H2')
    oh_idx=gas.species_index('OH')
    profiles=np.column_stack((f.grid, f.T, f.velocity, f.X[h2_idx,:], f.X[oh_idx,:]))
    np.savetxt(out/'flame_profiles.csv',profiles,delimiter=',',header='x_m,T_K,u_m_s,X_H2,X_OH',comments='')
    summary={'status':'Cantera FreeFlame chemical H2/O2/N2; no plasma',
             'mechanism':mechanism,'transport':'mixture-averaged',
             'flame_speed_m_s':speed,'peak_temperature_K':Tmax,'grid_points':len(f.grid),
             'width_m':width_m,'equivalence_ratio':phi,'pressure_atm':pressure_atm,
             'refine_slope':refine_slope,'refine_curve':refine_curve,
             'no_experimental_validation':True,
             'note':'Single mesh with automated refinement; requires independent grid and literature comparison.'}
    (out/'flame_summary.json').write_text(json.dumps(summary,indent=2))
    return summary


def run_flame_refinement_study(out, mechanism='h2o2.yaml',
                               tolerance=0.03, **kwargs):
    """Two independently solved Cantera flame grids for *numerical* sensitivity.

    Passing grid sensitivity is not agreement with a physical experiment.
    This requires Cantera, and is not executed in the offline release.
    """
    if not np.isfinite(tolerance) or not 0<tolerance<1:
        raise ValueError('Grid sensitivity tolerance must be in (0,1)')
    base=Path(out);base.mkdir(parents=True,exist_ok=True)
    coarse=run_flame(base/'coarse',mechanism=mechanism,
                     refine_slope=.12,refine_curve=.2,**kwargs)
    fine=run_flame(base/'fine',mechanism=mechanism,
                   refine_slope=.045,refine_curve=.08,**kwargs)
    error=abs(fine['flame_speed_m_s']-coarse['flame_speed_m_s'])/fine['flame_speed_m_s']
    report={'status':'GRID_SENSITIVITY_PASS' if error<=tolerance else 'GRID_SENSITIVITY_FAIL',
            'mechanism':mechanism,'coarse':coarse,'fine':fine,
            'relative_flame_speed_change':error,'tolerance':tolerance,
            'scientific_status':'Only two-grid solver sensitivity; NOT experimental flame validation'}
    (base/'grid_refinement_summary.json').write_text(json.dumps(report,indent=2))
    return report

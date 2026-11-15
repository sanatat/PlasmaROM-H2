"""Optional real two-term electron Boltzmann solve via external BOLOS.

Requires *user-supplied, properly sourced* complete BOLSIG+ collision file;
no cross-section data are embedded or invented. Output f is BOLOS' native
function and is NOT relabeled as a normalized EEDF probability PDF.
"""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
import numpy as np

def run_bolos(collision_file, out, *, source_url, mixture, field_Td=100.,gas_T_K=300.,
              upper_energy_eV=60.,grid_cells=500):
    if not source_url or not source_url.startswith(('https://','http://','doi:')):
        raise ValueError('Provide an attributable source URL/DOI for collision file')
    if not 0<field_Td<=1500 or not 100<=gas_T_K<=3000 or grid_cells<80 or upper_energy_eV<10:
        raise ValueError('Invalid Boltzmann solver settings')
    if not mixture or any(not np.isfinite(float(v)) or float(v)<=0 for v in mixture.values()) or not np.isclose(sum(mixture.values()),1):
        raise ValueError('Mixture must sum to 1 with strictly positive fractions')
    collision_file=Path(collision_file)
    if not collision_file.is_file():
        raise FileNotFoundError(collision_file)
    try:
        from bolos import parser, solver, grid
        from scipy.constants import k,e
    except ImportError as exc:
        raise RuntimeError('BOLOS unavailable; install separately with pip install bolos and compatible NumPy/SciPy') from exc
    gr=grid.LinearGrid(0,float(upper_energy_eV),grid_cells)
    b=solver.BoltzmannSolver(gr)
    with collision_file.open(encoding='utf-8') as fp:
        b.load_collisions(parser.parse(fp))
    if set(mixture)!=set(b.target):
        raise ValueError(f'Collision targets {sorted(b.target)} differ from mixture {sorted(mixture)}')
    for name,fraction in mixture.items():
        b.target[name].density=fraction
    b.kT=gas_T_K*k/e
    b.EN=field_Td*solver.TOWNSEND
    b.init()
    f=b.converge(b.maxwell(2.0),maxn=200,rtol=1e-5)
    if not np.isfinite(f).all() or np.any(f< -1e-12):
        raise RuntimeError('Boltzmann EEDF failed positivity/finite check')
    mean=float(b.mean_energy(f))
    if not np.isfinite(mean) or mean<=0:
        raise RuntimeError('Boltzmann mean energy nonphysical')
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    np.savetxt(out/'bolos_native_distribution.csv',np.column_stack([b.grid.c,f]),delimiter=',',
               header='energy_eV,bolos_native_eedf_not_probability_density',comments='')
    report={'status':'BOLOS_TWO_TERM_NUMERICAL_EEDF_NOT_DISCHARGE_SOLVER',
            'cross_section_source':source_url,'cross_section_sha256':hashlib.sha256(collision_file.read_bytes()).hexdigest(),
            'field_Td':field_Td,'temperature_K':gas_T_K,'mixture':mixture,
            'mean_energy_eV':mean,'mobility_times_neutral_density':float(b.mobility(f)),
            'diffusion_times_neutral_density':float(b.diffusion(f)),
            'note':'Native BOLOS EEDF, not probability_per_eV; no experimental validation and no plasma discharge.'}
    (out/'bolos_summary.json').write_text(json.dumps(report,indent=2))
    return report

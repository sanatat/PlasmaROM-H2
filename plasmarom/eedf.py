"""EEDF import and rate convolution with externally obtained electron distributions.

Input energy_eV,probability_per_eV is an energy-space *probability density*
normalized to 1. It is NOT BOLSIG+'s EEPF without the required conversion.
No Boltzmann solver or electric field calculation is included.
"""
from pathlib import Path
import csv
import numpy as np
from scipy.integrate import trapezoid
from scipy.constants import e as E_CHARGE, m_e


def validate_eedf(table, *, normalisation_tol=0.005):
    a=np.asarray(table,float)
    if (a.ndim != 2 or a.shape[0] < 5 or a.shape[1] != 2 or
        not np.all(np.isfinite(a)) or np.any(a < 0) or a[0,0] != 0 or
        not np.all(np.diff(a[:,0])>0)):
        raise ValueError('EEDF must be an increasing 0-based [energy_eV, probability_per_eV] table')
    integral=float(trapezoid(a[:,1],a[:,0]))
    if not 1-normalisation_tol<=integral<=1+normalisation_tol:
        raise ValueError(f'EEDF must integrate to unity (observed integral={integral:.6g})')
    return integral


def read_eedf_csv(path):
    with Path(path).open(newline='',encoding='utf-8') as f:
        rows=list(csv.DictReader(f))
    try:
        a=np.array([[float(z['energy_eV']),float(z['probability_per_eV'])] for z in rows])
    except (ValueError,TypeError,KeyError) as exc:
        raise ValueError('EEDF CSV must contain energy_eV,probability_per_eV') from exc
    validate_eedf(a)
    return a


def rate_from_eedf(cross_section, eedf, threshold_eV, *, max_outside_probability=1e-4):
    """Numerically evaluate <sigma v> [m³/s] from normalized energy PDF.

    Reject an EEDF whose probability above the known cross-section domain is
    too large. No extrapolation is allowed, and missing energy coverage is not
    silently treated as zero cross-section.
    """
    sigma=np.asarray(cross_section,float)
    p=np.asarray(eedf,float)
    integral = validate_eedf(p)
    if (sigma.ndim!=2 or sigma.shape[0]<3 or sigma.shape[1]!=2 or
        not np.all(np.isfinite(sigma)) or np.any(sigma<0) or
        sigma[0,0]!=0 or not np.all(np.diff(sigma[:,0])>0) or
        not np.isfinite(threshold_eV) or threshold_eV<0 or threshold_eV>=sigma[-1,0]):
        raise ValueError('Invalid cross-section table or threshold')
    if not 0<=max_outside_probability<1:
        raise ValueError('Invalid probability-tail limit')
    # Measure population in energy regions where cross sections are unknown.
    overlap_end=float(sigma[-1,0])
    grid=np.unique(np.r_[p[:,0],overlap_end])
    pdf=np.interp(grid,p[:,0],p[:,1],left=0,right=0)
    omitted=float(trapezoid(pdf[grid>=overlap_end],grid[grid>=overlap_end]))
    if omitted>max_outside_probability:
        raise ValueError(f'EEDF probability outside cross-section domain: {omitted:.5g}')
    valid=(p[:,0]<=overlap_end)
    energy=p[valid,0]
    pdf=p[valid,1]/integral
    if energy[-1]<overlap_end and overlap_end<=p[-1,0]:
        energy=np.r_[energy,overlap_end]
        pdf=np.r_[pdf,np.interp(overlap_end,p[:,0],p[:,1])/integral]
    sig=np.interp(energy,sigma[:,0],sigma[:,1])
    sig[energy<threshold_eV]=0
    if threshold_eV>0:
        energy=np.unique(np.r_[energy,threshold_eV])
        pdf=np.interp(energy,p[:,0],p[:,1])/integral
        sig=np.interp(energy,sigma[:,0],sigma[:,1])
        sig[energy<=threshold_eV]=0
    velocity=np.sqrt(2*E_CHARGE*energy/m_e)
    return float(trapezoid(sig*velocity*pdf,energy))


def maxwellian_pdf(energy_eV, te_eV):
    """Normalized Maxwellian ENERGY distribution [1/eV] on 0..infinity."""
    x=np.asarray(energy_eV,float)
    if (not np.isfinite(te_eV) or te_eV<=0 or not np.all(np.isfinite(x)) or np.any(x<0)):
        raise ValueError('Invalid energy or Te')
    return 2/np.sqrt(np.pi)*np.sqrt(x)*np.exp(-x/te_eV)/te_eV**1.5


def run_eedf_rate(out, eedf_csv, cross_section_csv, source_eedf, source_sigma, threshold_eV):
    """Read traceable external inputs and save a single electron-impact rate.

    Provenance strings must give a dataset URL, accession, or citation; their
    existence is not independently checked by this software.
    """
    import hashlib
    import json
    from .plasma_states import read_cross_section_csv
    if not source_eedf or not source_sigma:
        raise ValueError('Both data sources must be documented')
    eedf=read_eedf_csv(eedf_csv)
    sigma=read_cross_section_csv(cross_section_csv)
    rate=rate_from_eedf(sigma,eedf,threshold_eV)
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    def hash_file(f):return hashlib.sha256(Path(f).read_bytes()).hexdigest()
    result={'rate_coefficient_m3_s':rate,'threshold_eV':float(threshold_eV),
            'eedf_integral':validate_eedf(eedf),
            'source_eedf':source_eedf,'source_cross_section':source_sigma,
            'eedf_sha256':hash_file(eedf_csv),'cross_section_sha256':hash_file(cross_section_csv),
            'status':'user-supplied externally computed EEDF; no Boltzmann or discharge simulation'}
    (out/'external_eedf_rate.json').write_text(json.dumps(result,indent=2))
    return result

"""Interpolate externally solved EEDF tables across reduced field E/N [Td].

This module does NOT solve the Boltzmann equation or predict E/N from a
self-consistent discharge. Only applies to a documented external EEDF grid;
interpolation is a numerical convenience with unquantified model error.
"""
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from .eedf import read_eedf_csv, validate_eedf, rate_from_eedf
from .plasma_states import read_cross_section_csv


def read_field_manifest(path):
    """CSV: EoverN_Td,eedf_csv,source. Paths relative to manifest folder."""
    manifest=Path(path)
    with manifest.open(newline='',encoding='utf-8') as fh:
        rows=list(csv.DictReader(fh))
    if len(rows)<2:
        raise ValueError('Need at least two independently supplied E/N EEDF states')
    data=[]
    for row in rows:
        try:
            field=float(row['EoverN_Td'])
            name=row['eedf_csv']
            source=row['source'].strip()
        except (ValueError,KeyError,TypeError,AttributeError) as ex:
            raise ValueError('Manifest needs EoverN_Td,eedf_csv,source') from ex
        if not np.isfinite(field) or field<=0 or not name or not source:
            raise ValueError('Positive E/N, EEDF file, and provenance are required')
        file_path=(manifest.parent/name).resolve()
        if not file_path.is_file():
            raise FileNotFoundError(file_path)
        pdf=read_eedf_csv(file_path)
        data.append((field,pdf,source,hashlib.sha256(file_path.read_bytes()).hexdigest()))
    data.sort(key=lambda d:d[0])
    if len(set(z[0] for z in data))!=len(data):
        raise ValueError('Reduced field E/N entries must be unique')
    reference_E=data[0][1][:,0]
    if any(not np.array_equal(z[1][:,0],reference_E) for z in data[1:]):
        raise ValueError('For auditable interpolation, all EEDF grids must match exactly')
    return data


def interpolate_field_eedf(entries, EoverN_Td):
    """Piecewise interpolation of probability densities in logarithmic E/N."""
    fields=np.array([x[0] for x in entries],float)
    value=float(EoverN_Td)
    if not np.isfinite(value) or value<fields[0] or value>fields[-1]:
        raise ValueError('Requested E/N outside measured external EEDF grid; extrapolation forbidden')
    hi=int(np.searchsorted(fields,value,side='left'))
    if hi==0 or fields[hi]==value:
        pdf=entries[hi][1].copy()
        provenance=[entries[hi][2]]
    else:
        lo=hi-1
        a=(np.log(value)-np.log(fields[lo]))/(np.log(fields[hi])-np.log(fields[lo]))
        pdf=entries[lo][1].copy()
        pdf[:,1]=(1-a)*entries[lo][1][:,1]+a*entries[hi][1][:,1]
        provenance=[entries[lo][2],entries[hi][2]]
    validate_eedf(pdf)
    return pdf,provenance


def run_field_rate(out, manifest, sigma_csv, sigma_source, threshold_eV, field_Td):
    if not sigma_source or not sigma_source.strip():
        raise ValueError('Cross-section provenance required')
    entries=read_field_manifest(manifest)
    eedf,source=interpolate_field_eedf(entries,field_Td)
    sigma=read_cross_section_csv(sigma_csv)
    rate=rate_from_eedf(sigma,eedf,threshold_eV)
    res={'EoverN_Td':float(field_Td),'rate_coefficient_m3_s':rate,
         'source_eedf_fields':source,'source_cross_section':sigma_source,
         'all_field_values_Td':[x[0] for x in entries],
         'eedf_sha256_by_field':{str(x[0]):x[3] for x in entries},
         'cross_section_sha256':hashlib.sha256(Path(sigma_csv).read_bytes()).hexdigest(),
         'scientific_status':'Interpolated externally supplied EEDF only; NOT solved Boltzmann or validated PAC',
         'interpolation':'log(E/N) linear blend of probability densities, no extrapolation'}
    output=Path(out);output.mkdir(parents=True,exist_ok=True)
    (output/'field_dependent_rate.json').write_text(json.dumps(res,indent=2),encoding='utf-8')
    return res

"""Compare imported independent measurements with numerical predictions.

Never treat a synthetic table as experimental validation. Every measurement
needs an explicit source, units, uncertainty and independent set label.
"""
import csv
import json
from pathlib import Path
import numpy as np

def evaluate_experiment(measurements, predictions, output):
    def load(f):
        with Path(f).open(encoding='utf-8',newline='') as fh:
            return list(csv.DictReader(fh))
    ms=load(measurements);ps=load(predictions)
    if not ms or not ps:raise ValueError('Empty comparison input')
    index={}
    for row in ps:
        key=(row.get('case_id'),row.get('observable'))
        if not all(key) or key in index:raise ValueError('Duplicate/missing prediction key')
        if not np.isfinite(float(row['predicted_value'])):raise ValueError('Nonfinite prediction')
        index[key]=row
    errors=[];sources=set();observed_keys=set()
    for m in ms:
        key=(m.get('case_id'),m.get('observable'))
        if not all(key):
            raise ValueError('Missing measurement case_id or observable')
        if key in observed_keys:
            raise ValueError('Duplicate measurement (case_id, observable). Assign each genuinely independent replicate a distinct case_id and provide a corresponding prediction; do not silently double-count it.')
        observed_keys.add(key)
        if key not in index:raise ValueError('Missing model prediction for observed case')
        if m.get('data_origin')!='MEASURED_INDEPENDENT' or not m.get('source_doi_or_url') or not m.get('units'):
            raise ValueError('Measurement provenance incomplete or not declared independent')
        v=float(m['measured_value']);u=float(m['uncertainty_1sigma']);p=float(index[key]['predicted_value'])
        if not np.isfinite([v,u,p]).all() or u<=0 or m['units']!=index[key]['units']:
            raise ValueError('Uncertainty/units/values invalid')
        errors.append({'case_id':key[0],'observable':key[1], 'z_score':(p-v)/u,
                       'abs_error':abs(p-v),'units':m['units'],'source':m['source_doi_or_url']})
        sources.add(m['source_doi_or_url'])
    z=np.array([e['z_score'] for e in errors])
    report={'status':'COMPARISON_TO_USER_PROVIDED_INDEPENDENT_MEASUREMENTS_NOT_AUTOMATICALLY_CERTIFIED',
            'n_measurements':len(errors),'rmse_in_uncertainty_units':float(np.sqrt(np.mean(z*z))),
            'max_abs_z_score':float(np.max(np.abs(z))),'sources':sorted(sources),'cases':errors,
            'note':'Source URL is recorded but neither study credibility nor sample independence is verified by software.'}
    out=Path(output);out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps(report,indent=2))
    return report

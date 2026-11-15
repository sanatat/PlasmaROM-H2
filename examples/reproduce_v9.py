"""Standalone, SYNTHETIC-only numerical and scenario-protocol demonstration."""
from pathlib import Path
import json
import platform
import sys
# Allow an uninstalled checkout to reproduce numerical verification.
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from plasmarom.acceptance import _conditions,_overlap
from plasmarom.verification import convergence_study


def main():
    root=Path(__file__).resolve().parents[1]
    scenarios=json.loads((root/'examples/cases/h2_composition_sweep.json').read_text())
    train=_conditions(scenarios['training'])
    holdout=_conditions(scenarios['holdout'])
    if train is None or holdout is None or _overlap(train,holdout):
        raise RuntimeError('Scenario protocol invalid or overlapping')
    out=root/'results_v9'
    out.mkdir(parents=True,exist_ok=True)
    numeric=convergence_study(out=out)
    report={
        'software_version':'0.9.0',
        'python':sys.version.split()[0],
        'platform':platform.platform(),
        'training_cases':len(train),
        'heldout_cases':len(holdout),
        'composition_aware_split':'PASSED',
        'numerical_convergence':{k:{'observed_order':float(v['observed_order'])} for k,v in numeric.items()},
        'physical_plasma_validation':False,
        'real_cantera_kinetics_verified_by_this_script':False,
        'validated_skeletal_mechanism':False,
        'notice':'Only synthetic and manufactured numerical tests; Cantera requires separate run.'
    }
    dest=out/'REPRODUCE_V9.json'
    dest.write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps(report,indent=2))
    return report

if __name__=='__main__':main()

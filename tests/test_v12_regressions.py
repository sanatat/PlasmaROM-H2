"""Regression tests for four concrete defects identified in v11.

All synthetic evidence is unit-test input, not real PAC validation.
"""
import csv
from pathlib import Path
import numpy as np
import pytest

from plasmarom.ammonia_h2 import check_initial_mixture
from plasmarom.experiment_compare import evaluate_experiment
from plasmarom.energy_graph import energy_aware_screen


def test_co_counts_in_stoichiometric_oxygen():
    mix = {'NH3': 1, 'H2': 1, 'CO': 1, 'O2': 1.25, 'N2': 4.7}
    assert check_initial_mixture(mix) == pytest.approx(1.4)


def test_previous_co_free_case_unchanged():
    assert check_initial_mixture({'NH3':1,'H2':1,'O2':1.25,'N2':4.7}) == pytest.approx(1)


def test_unknown_fuel_fails_instead_of_silent_wrong_phi():
    with pytest.raises(ValueError, match='unsupported species'):
        check_initial_mixture({'NH3':1,'H2':1,'CH4':1,'O2':1.25,'N2':4.7})


def test_negative_or_nonfinite_fuel_rejected():
    for value in (-1., float('nan'), float('inf'), True):
        with pytest.raises(ValueError):
            check_initial_mixture({'NH3':1,'H2':1,'CO':value,'O2':1.25,'N2':4.7})


def write_csv(path, columns, rows):
    with path.open('w', newline='', encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=columns)
        w.writeheader()
        w.writerows(rows)


def experiment_rows(tmp_path, cases):
    m=tmp_path/'meas.csv';p=tmp_path/'pred.csv'
    measurements=[dict(case_id=c,observable=o,measured_value='2',uncertainty_1sigma='0.1',units='s',source_doi_or_url='https://doi.org/10.1000/test',data_origin='MEASURED_INDEPENDENT') for c,o in cases]
    predictions=[dict(case_id=c,observable=o,predicted_value='2.1',units='s') for c,o in set(cases)]
    write_csv(m,['case_id','observable','measured_value','uncertainty_1sigma','units','source_doi_or_url','data_origin'],measurements)
    write_csv(p,['case_id','observable','predicted_value','units'],predictions)
    return m,p


def test_duplicate_measurement_is_rejected(tmp_path):
    m,p=experiment_rows(tmp_path,[('c1','ignition'),('c1','ignition')])
    with pytest.raises(ValueError,match='Duplicate measurement'):
        evaluate_experiment(m,p,tmp_path/'result.json')
    assert not (tmp_path/'result.json').exists()


def test_different_observables_same_case_are_allowed(tmp_path):
    m,p=experiment_rows(tmp_path,[('c1','ignition'),('c1','flame')])
    report=evaluate_experiment(m,p,tmp_path/'result.json')
    assert report['n_measurements']==2


def test_explicit_replicates_need_distinct_case_ids(tmp_path):
    m,p=experiment_rows(tmp_path,[('c1_rep1','ignition'),('c1_rep2','ignition')])
    assert evaluate_experiment(m,p,tmp_path/'result.json')['n_measurements']==2


def test_superelastic_gain_is_supported_and_retained():
    names=['target','excitation','superelastic']
    g=np.zeros((3,3))
    rates=np.tile([1.,100.,50.],(11,1))
    d=energy_aware_screen(g,['target'],rates,[0.,2.,-3.],np.linspace(0,1,11),names,graph_threshold=0.9,min_energy_fraction=0.99)
    assert set(d['selected'])==set(names)
    assert d['retained_loss_fraction']>=0.99
    assert d['retained_gain_fraction']>=0.99
    assert d['integrated_electron_energy_eV_units']['superelastic']<0
    assert d['net_energy_exchange_eV_units']==pytest.approx(50.)


def test_superelastic_only_no_false_zero_total():
    d=energy_aware_screen(np.zeros((2,2)),['target'],np.ones((3,2)),[0.,-5.],[0,.5,1],['target','gain'],min_energy_fraction=.95)
    assert 'gain' in d['selected']
    assert d['retained_energy_fraction']==pytest.approx(1.)
    assert d['retained_loss_fraction'] is None
    assert d['retained_gain_fraction']==pytest.approx(1.)


def test_energy_graph_rejects_negative_channel_rates():
    with pytest.raises(ValueError):
        energy_aware_screen(np.zeros((2,2)),['a'],np.array([[1,1],[1,-1.]]),[-1.,2.],[0.,1.],['a','b'])


def test_workflow_covers_changes_to_all_model_files():
    config=Path(__file__).parents[1]/'.github/workflows/v12-real-chemistry.yml'
    workflow=config.read_text()
    assert 'pull_request:' in workflow
    assert "'plasmarom/**'" in workflow
    assert "'tests/**'" in workflow
    assert "'examples/**'" in workflow
    assert "'pyproject.toml'" in workflow
    assert 'test_v12_regressions.py' in workflow
    assert not (config.parent/'v11-real-chemistry.yml').exists()

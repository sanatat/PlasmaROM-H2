"""Version 9 regression tests: complete provenance, independent replay, mixtures.

Fake Cantera metadata and replay objects test fail-closed control FLOW only;
actual reaction-kinetics integration needs the optional cantera package in CI.
"""
from dataclasses import dataclass
from types import SimpleNamespace
import pytest
from plasmarom import acceptance
from plasmarom.cantera_demo import run_skeletal
from test_v7_evidence_integrity import sample_summary, fake_verified_metadata


def _patch_verified(monkeypatch, summary):
    monkeypatch.setattr(acceptance,'_load_mechanism_metadata',fake_verified_metadata)
    monkeypatch.setattr(acceptance,'_replay_holdout_metrics',
                        lambda *args:[{k:v for k,v in row.items() if k!='condition_T_K_P_atm'}
                                      for row in summary['holdout_metrics']])
    monkeypatch.setattr(acceptance,'_replay_training_observations',
                        lambda *args:[{k:v for k,v in row.items() if k!='condition_T_K_P_atm'}
                                      for row in summary['training_ignition_observations']])


def test_valid_policy_fixture_still_passes_with_all_stubs(tmp_path,monkeypatch):
    data=sample_summary(tmp_path)
    _patch_verified(monkeypatch,data)
    assert acceptance.assess_skeletal_summary(data,tmp_path)['accepted'] is True


def test_v9_forged_training_rejected_even_with_self_consistent_counts(tmp_path,monkeypatch):
    data=sample_summary(tmp_path)
    real_records=[dict(x) for x in data['training_ignition_observations']]
    _patch_verified(monkeypatch,data)
    monkeypatch.setattr(acceptance,'_replay_training_observations',
                        lambda *args:[{k:v for k,v in x.items() if k!='condition_T_K_P_atm'}
                                      for x in real_records])
    # Keep the aggregate count plausible while falsifying *which* training case ignites.
    data['training_ignition_observations'][0]['reference_ignition_observed']=True
    data['training_ignition_observations'][1]['reference_ignition_observed']=False
    result=acceptance.assess_skeletal_summary(data,tmp_path)
    assert not result['accepted']
    assert 'independent_training_replay_mismatch' in result['reasons']


def test_v9_training_replay_failure_fails_closed(tmp_path,monkeypatch):
    data=sample_summary(tmp_path)
    _patch_verified(monkeypatch,data)
    def fail(*args): raise RuntimeError('reactor integration failed')
    monkeypatch.setattr(acceptance,'_replay_training_observations',fail)
    result=acceptance.assess_skeletal_summary(data,tmp_path)
    assert not result['accepted']
    assert 'independent_training_replay_failed' in result['reasons']


def test_v9_missing_training_clocks_rejected(tmp_path,monkeypatch):
    data=sample_summary(tmp_path)
    _patch_verified(monkeypatch,data)
    del data['training_ignition_observations'][0]['max_time_step_s']
    result=acceptance.assess_skeletal_summary(data,tmp_path)
    assert not result['accepted']
    assert 'independent_training_replay_mismatch' in result['reasons']


def test_v9_modified_thermodynamics_rejected(tmp_path,monkeypatch):
    data=sample_summary(tmp_path)
    original=fake_verified_metadata(None,'h2o2.yaml')
    original['reduced_species_fingerprints']['X1']='NASA7:changed'
    monkeypatch.setattr(acceptance,'_load_mechanism_metadata',lambda *args:original)
    result=acceptance.assess_skeletal_summary(data,tmp_path)
    assert not result['accepted']
    assert 'modified_species_thermo_or_transport' in result['reasons']


def test_v9_no_species_fingerprint_fails_closed(tmp_path,monkeypatch):
    data=sample_summary(tmp_path)
    metadata=fake_verified_metadata(None,'h2o2.yaml')
    del metadata['source_species_fingerprints']
    monkeypatch.setattr(acceptance,'_load_mechanism_metadata',lambda *args:metadata)
    result=acceptance.assess_skeletal_summary(data,tmp_path)
    assert 'missing_full_species_identity_evidence' in result['reasons']


def test_v9_species_thermo_parameters_are_in_fingerprint():
    def sp(coeff):
        return SimpleNamespace(name='H2O',composition={'H':2,'O':1},charge=0,size=1,
            input_data={'name':'H2O','thermo':{'model':'NASA7','data':[[1,coeff]]}},
            thermo=SimpleNamespace(coeffs=[1.,coeff],min_temp=200.,max_temp=3500.))
    a=acceptance._species_fingerprint(sp(2.0))
    b=acceptance._species_fingerprint(sp(2.1))
    assert a!=b


def test_v9_duplicate_mixtures_with_scaled_coefficients_rejected():
    cases=[[1100.,1.,'H2:2,O2:1,N2:4'],[1100.,1.,'N2:8, H2:4, O2:2']]
    assert acceptance._conditions(cases) is None


def test_v9_mixed_composition_cases_are_distinct():
    a=acceptance._conditions([[1100.,1.,'H2:2,O2:1,N2:4']])
    b=acceptance._conditions([[1100.,1.,'H2:1,O2:1,N2:4']])
    assert a is not None and b is not None
    assert not acceptance._overlap(a,b)


@pytest.mark.parametrize('composition', ['H2:-1,O2:1','H2:nan,O2:1',
                                        'H2:1,H2:2','H2:1,,O2:2','',
                                        'H2:1,O2:0','H2:inf'])
def test_v9_bad_mixtures_rejected_before_cantera(tmp_path,composition):
    with pytest.raises(ValueError,match='disjoint and unique'):
        run_skeletal(tmp_path,training=((1100.,1.,composition),),holdout=((1200.,1.),))
    assert not list(tmp_path.iterdir())


def test_v9_mixture_protocol_must_match_summary(tmp_path,monkeypatch):
    data=sample_summary(tmp_path)
    _patch_verified(monkeypatch,data)
    data['holdout_conditions'][0].append('H2:1,O2:1,N2:4')
    data['holdout_metrics'][0]['condition_T_K_P_atm'].append('H2:1,O2:1,N2:4')
    result=acceptance.assess_skeletal_summary(data,tmp_path)
    assert 'composition_protocol_mismatch' in result['reasons']


def test_v9_distinct_mixture_cases_are_validated_prior_to_optional_import():
    a=acceptance._conditions([[1100.,1.,'H2:2,O2:1,N2:4']])
    b=acceptance._conditions([[1100.,1.,'H2:3,O2:1,N2:4']])
    assert len(a)==1 and len(b)==1
    assert not acceptance._overlap(a,b)


def test_v9_mixture_species_provenance_validation():
    assert acceptance._composition_signature(' NH3:1, H2:3, O2:2, N2:8 ') is not None
    assert acceptance._composition_signature(' NH3:1, H2:3, O2:2, N2:8 ') == \
           acceptance._composition_signature('N2:16,O2:4,H2:6,NH3:2')

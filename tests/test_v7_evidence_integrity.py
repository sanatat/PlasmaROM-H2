"""Falsification tests for acceptance, data isolation and artifact verification.

Metadata stubs here prove only the policy logic; real integration tests are
optional and require importing Cantera.
"""
import copy
import hashlib
import json

import pytest

import plasmarom.acceptance as acceptance
from plasmarom.cantera_demo import run_skeletal
from plasmarom.channel_reduction import evaluate_channel_pruning


def sample_summary(tmp_path):
    p = tmp_path / 'holdout_pass_skeletal.yaml'
    p.write_text('INTENTIONALLY NOT A REAL CANTERA MECHANISM', encoding='utf-8')
    return {
        'status':'HOLDOUT_PASS',
        'full_species':10, 'reduced_species':8,
        'full_reactions':30, 'reduced_reactions':20,
        'source_mechanism':'h2o2.yaml',
        'source_mechanism_sha256':'test-source-hash-fixture',
        'reactor_t_end_s':.006, 'reactor_n_snapshots':80,
        'initial_composition':'H2:2,O2:1,N2:4',
        'ignition_definition':'peak_dTdt_if_threshold_resolved_with_20K_rise',
        'training_conditions':[[1050.,1.],[1200.,1.]],
        'training_ignition_observations':[
            {'condition_T_K_P_atm':[1050.,1.],'reference_ignition_observed':False,
             'threshold_crossing_time_s':None,'peak_dTdt_time_s':None,
             'max_time_step_s':.006/79},
            {'condition_T_K_P_atm':[1200.,1.],'reference_ignition_observed':True,
             'threshold_crossing_time_s':.002,'peak_dTdt_time_s':.0018,
             'max_time_step_s':.006/79}],
        'ignition_training_cases_observed':1,
        'holdout_conditions':[[1150.,1.]],
        'holdout_metrics':[
            {'condition_T_K_P_atm':[1150.,1.],
             'reference_ignition_observed':True,
             'candidate_ignition_observed':True,
             'ignition_classification_matches':True,
             'ignition_error_evaluable':True,
             'ignition_proxy_rel_error':.02,
             'ignition_time_definition':'peak_dTdt_if_threshold_resolved',
             'temperature_threshold_rise_K':20.,
             'reference_threshold_crossing_time_s':.002,
             'candidate_threshold_crossing_time_s':.0021,
             'reference_peak_dTdt_time_s':.0018,
             'candidate_peak_dTdt_time_s':.00185,
             'reference_max_time_step_s':.006/79,
             'candidate_max_time_step_s':.006/79,
             'T_K_normalized_max_abs_error':.03,
             'OH_Y_normalized_max_abs_error':.01,
             'H2O_Y_normalized_max_abs_error':.01}],
        'ignition_holdout_cases_observed':1,
        'holdout_ignition_classifications_match':True,
        'holdout_max_normalized_error':.03,
        'holdout_tolerance':.12,
        'accepted_with_ignition_evidence':True,
        'artifact':'holdout_pass_skeletal.yaml',
        'artifact_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),
    }


def fake_verified_metadata(path,source):
    assert source == 'h2o2.yaml'
    return {'full_species':10, 'reduced_species':8,
            'full_reactions':30,'reduced_reactions':20,
            'source_species_names':tuple(f'X{k}' for k in range(10)),
            'reduced_species_names':tuple(f'X{k}' for k in range(8)),
            'source_species_composition':{f'X{k}':{'H':float(k+1)} for k in range(10)},
            'reduced_species_composition':{f'X{k}':{'H':float(k+1)} for k in range(8)},
            'source_species_fingerprints':{f'X{k}':f'meta-{k}' for k in range(10)},
            'reduced_species_fingerprints':{f'X{k}':f'meta-{k}' for k in range(8)},
            'source_sha256':'test-source-hash-fixture',
            'source_reaction_signatures':tuple(f'reaction-{j}' for j in range(30)),
            'reduced_reaction_signatures':tuple(f'reaction-{j}' for j in range(20))}


def test_fake_yaml_does_not_pass_even_with_perfect_self_report(tmp_path):
    summary = sample_summary(tmp_path)
    check = acceptance.assess_skeletal_summary(summary,tmp_path)
    assert check['accepted'] is False
    assert 'invalid_or_unreadable_cantera_mechanism' in check['reasons'] or 'cantera_unavailable_cannot_validate_mechanism' in check['reasons']


def test_consistent_evidence_requires_a_verified_phase(tmp_path,monkeypatch):
    summary = sample_summary(tmp_path)
    monkeypatch.setattr(acceptance,'_load_mechanism_metadata',fake_verified_metadata)
    monkeypatch.setattr(acceptance,'_replay_holdout_metrics',lambda *args: summary['holdout_metrics'])
    monkeypatch.setattr(acceptance,'_replay_training_observations',
                        lambda *args:[{k:v for k,v in row.items() if k!='condition_T_K_P_atm'}
                                      for row in summary['training_ignition_observations']])
    assert acceptance.assess_skeletal_summary(summary,tmp_path)['accepted'] is True
    assert acceptance.assess_skeletal_summary(summary,tmp_path,require_artifact=False)['accepted'] is False


@pytest.mark.parametrize('mutate,reason',[
    (lambda s:s['holdout_metrics'][0].update(candidate_ignition_observed=False),
     'contradictory_ignition_classification'),
    (lambda s:s['holdout_metrics'][0].update(ignition_proxy_rel_error=-.1),
     'invalid_ignition_error'),
    (lambda s:s.update(ignition_holdout_cases_observed=99),
     'inconsistent_holdout_ignition_count'),
    (lambda s:s.update(ignition_training_cases_observed=99),
     'inconsistent_training_ignition_count'),
    (lambda s:s.update(holdout_max_normalized_error=.001),
     'aggregate_error_inconsistent_with_cases'),
    (lambda s:s.update(holdout_conditions=[[1050.,1.]]),
     'training_holdout_overlap'),
    (lambda s:s['holdout_metrics'][0].update(condition_T_K_P_atm=[1300,1]),
     'holdout_condition_record_mismatch'),
    (lambda s:s['training_ignition_observations'][0].update(reference_ignition_observed=True),
     'inconsistent_training_ignition_count'),
    (lambda s:s.update(holdout_ignition_classifications_match=False),
     'inconsistent_aggregate_ignition_status'),
    (lambda s:s.update(accepted_with_ignition_evidence=False),
     'inconsistent_claimed_acceptance'),
    (lambda s:s.update(reduced_reactions=30),
     'no_nontrivial_mechanism_reduction'),
    (lambda s:s.update(training_conditions=[[1050.,1.],[1050.,1.]]),
     'invalid_training_conditions'),
])
def test_forged_or_inconsistent_summary_rejected(tmp_path,monkeypatch,mutate,reason):
    s=sample_summary(tmp_path)
    monkeypatch.setattr(acceptance,'_load_mechanism_metadata',fake_verified_metadata)
    mutate(s)
    out=acceptance.assess_skeletal_summary(s,tmp_path)
    assert not out['accepted'],out
    assert reason in out['reasons'],out


def test_mechanism_metadata_must_match_declared_counts_and_species(tmp_path,monkeypatch):
    sample=sample_summary(tmp_path)
    inconsistent=fake_verified_metadata(None,'h2o2.yaml')
    inconsistent['reduced_reactions']=19
    monkeypatch.setattr(acceptance,'_load_mechanism_metadata',lambda *args:inconsistent)
    out=acceptance.assess_skeletal_summary(sample,tmp_path)
    assert 'cantera_reduced_reactions_count_mismatch' in out['reasons']
    inconsistent['reduced_reactions']=20
    inconsistent['reduced_species_names']=('NOT_IN_PARENT',)
    out=acceptance.assess_skeletal_summary(sample,tmp_path)
    assert 'reduced_species_not_subset_of_source' in out['reasons']


def test_corrupted_artifact_hash_is_rejected(tmp_path,monkeypatch):
    s=sample_summary(tmp_path)
    monkeypatch.setattr(acceptance,'_load_mechanism_metadata',fake_verified_metadata)
    (tmp_path/'holdout_pass_skeletal.yaml').write_text('changed')
    assert 'skeletal_yaml_hash_mismatch' in acceptance.assess_skeletal_summary(s,tmp_path)['reasons']


@pytest.mark.parametrize('training,holdout',[
    ((500.,),(500.,)), ((500.,),(500.00000000001,)),
    ((500.,500.),(700.,)), ((500.,),(700.,700.))])
def test_energy_channel_train_holdout_leakage_fails_before_simulation(tmp_path,training,holdout):
    with pytest.raises(ValueError,match='overlap|Duplicate'):
        evaluate_channel_pruning(tmp_path,training_powers=training,holdout_powers=holdout)
    assert not (tmp_path/'channel_reduction_summary.json').exists()


@pytest.mark.parametrize('training,holdout',[
    (((1050.,1.),),((1050.,1.),)),
    (((1050.,1.),),((1050.00000000001,1.),)),
    (((1050.,1.),(1050.,1.)),((1150.,1.),)),
    (((1050.,1.),),((1150.,1.),(1150.,1.)))])
def test_cantera_train_holdout_leakage_fails_before_optional_import(tmp_path,training,holdout):
    with pytest.raises(ValueError,match='disjoint and unique'):
        run_skeletal(tmp_path,training=training,holdout=holdout)
    assert not list(tmp_path.iterdir())


def test_acceptance_file_roundtrip_failure_is_reported(tmp_path):
    report=sample_summary(tmp_path)
    file=tmp_path/'skeletal_summary.json'
    file.write_text(json.dumps(report))
    result=acceptance.assess_skeletal_file(file,tmp_path/'decision.json')
    assert not result['accepted']
    assert json.loads((tmp_path/'decision.json').read_text())==result


def test_rejected_scientific_gate_cannot_leave_holdout_pass(tmp_path,monkeypatch):
    from plasmarom.cantera_demo import _finalize_scientific_acceptance
    from plasmarom import acceptance
    report=sample_summary(tmp_path)
    monkeypatch.setattr(acceptance,'assess_skeletal_summary',lambda *args,**kwargs:{'accepted':False,'reasons':['invalid_artifact']})
    _finalize_scientific_acceptance(report,tmp_path)
    assert report['status']=='FAILED_OR_INCONCLUSIVE_HOLDOUT'
    assert report['accepted_with_ignition_evidence'] is False
    assert report['acceptance_evidence']['reasons']==['invalid_artifact']
    assert report['artifact']=='candidate_skeletal.yaml'
    assert not (tmp_path/'holdout_pass_skeletal.yaml').exists()


def test_v8_modified_reaction_fingerprint_fails_even_when_counts_match(tmp_path,monkeypatch):
    summary=sample_summary(tmp_path)
    metadata=fake_verified_metadata(None,'h2o2.yaml')
    metadata['reduced_reaction_signatures']=('tampered-reaction',)+tuple(metadata['reduced_reaction_signatures'][1:])
    monkeypatch.setattr(acceptance,'_load_mechanism_metadata',lambda *args:metadata)
    monkeypatch.setattr(acceptance,'_replay_holdout_metrics',lambda *args:summary['holdout_metrics'])
    result=acceptance.assess_skeletal_summary(summary,tmp_path)
    assert not result['accepted']
    assert 'modified_or_foreign_reaction_parameters' in result['reasons']


def test_v8_source_mechanism_hash_is_required(tmp_path,monkeypatch):
    summary=sample_summary(tmp_path)
    monkeypatch.setattr(acceptance,'_load_mechanism_metadata',fake_verified_metadata)
    monkeypatch.setattr(acceptance,'_replay_holdout_metrics',lambda *args:summary['holdout_metrics'])
    summary['source_mechanism_sha256']='other-source-hash'
    result=acceptance.assess_skeletal_summary(summary,tmp_path)
    assert not result['accepted']
    assert 'source_mechanism_hash_mismatch' in result['reasons']


def test_v8_replay_rejects_self_consistent_but_forged_metrics(tmp_path,monkeypatch):
    summary=sample_summary(tmp_path)
    monkeypatch.setattr(acceptance,'_load_mechanism_metadata',fake_verified_metadata)
    replay=[dict(summary['holdout_metrics'][0])]
    replay[0]['T_K_normalized_max_abs_error']=0.07
    monkeypatch.setattr(acceptance,'_replay_holdout_metrics',lambda *args:replay)
    result=acceptance.assess_skeletal_summary(summary,tmp_path)
    assert not result['accepted']
    assert 'independent_reactor_replay_mismatch' in result['reasons']


def test_v8_replay_failure_closes_gate(tmp_path,monkeypatch):
    summary=sample_summary(tmp_path)
    monkeypatch.setattr(acceptance,'_load_mechanism_metadata',fake_verified_metadata)
    def broken(*args):
        raise RuntimeError('reactor divergence')
    monkeypatch.setattr(acceptance,'_replay_holdout_metrics',broken)
    result=acceptance.assess_skeletal_summary(summary,tmp_path)
    assert not result['accepted']
    assert 'independent_reactor_replay_failed' in result['reasons']

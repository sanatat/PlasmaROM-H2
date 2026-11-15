"""v10 regressions: collider-aware gross flux, trace OH, pinned plans and replay.

Cantera behaviour remains unverified here without Cantera. Fake metadata proves
only correct decision branches, not actual chemical-model validity.
"""
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import pytest
from plasmarom import acceptance
from plasmarom.protocol import (load_pinned_protocol, verify_protocol_parameters,
                                 verify_summary_protocol)
from plasmarom.cantera_demo import _reaction_participation_matrix,_curve_error
from test_v7_evidence_integrity import sample_summary, fake_verified_metadata


class FakePhase:
    species_names=('H2','O2','N2','AR')
    reactant_stoich_coeffs=np.asarray([[1.,0.],[1.,1.],[0.,0.],[0.,0.]])
    product_stoich_coeffs=np.asarray([[0.,0.],[0.,0.],[0.,1.],[0.,0.]])
    def __init__(self, third):
        self.third=third
    def reactions(self):
        return [SimpleNamespace(third_body=self.third,orders={}),
                SimpleNamespace(third_body=None,orders={})]
    def species_index(self,name):
        return self.species_names.index(name)


def test_colliders_in_gross_and_net_participation_are_consistent():
    third=SimpleNamespace(collider='M',efficiencies={'AR':3.},default_efficiency=0.)
    phase=FakePhase(third)
    participation=_reaction_participation_matrix(phase)
    assert bool(participation[3,0]) # AR is nonstoichiometric collider
    assert not bool(participation[3,1])
    assert bool(participation[0,0])
    assert bool(participation[1,0])


def test_collider_changes_gross_flux_graph_not_only_boolean_matrix():
    from plasmarom.reduction import gross_flux_interaction_coefficients
    phase=FakePhase(SimpleNamespace(collider='M',efficiencies={'AR':2.},
                                    default_efficiency=0.))
    nu=phase.product_stoich_coeffs-phase.reactant_stoich_coeffs
    p=_reaction_participation_matrix(phase)
    no_collider=p.copy(); no_collider[3,0]=False
    with_collider=gross_flux_interaction_coefficients(
        nu, np.array([4.,0.]),np.zeros(2),p)
    without=gross_flux_interaction_coefficients(
        nu, np.array([4.,0.]),np.zeros(2),no_collider)
    assert with_collider[0,3]>0.9
    assert without[0,3]==0


def test_generic_M_nonzero_default_includes_every_species():
    third=SimpleNamespace(collider='M',efficiencies={},default_efficiency=1.)
    phase=FakePhase(third)
    participation=_reaction_participation_matrix(phase)
    assert np.all(participation[:,0])
    assert not participation[3,1]


def test_explicit_third_body_is_included_even_without_stoichiometry():
    third=SimpleNamespace(collider='AR',efficiencies={},default_efficiency=0.)
    assert _reaction_participation_matrix(FakePhase(third))[3,0]


def test_trace_OH_metric_catches_large_relative_deviation():
    t=np.linspace(0,1,7)
    def curve(oh):
        return {'T_K':np.linspace(300,330,7),'OH_Y':np.full(7,oh),
                'H2O_Y':np.full(7,1e-5),'ignition_observed':False,
                'threshold_crossing_time_s':None,'peak_dTdt_time_s':None,
                'max_time_step_s':1/6,'temperature_threshold_rise_K':20.,
                'ignition_time_s':None}
    metrics=_curve_error(curve(1e-12),curve(5e-9))
    assert metrics['OH_Y_normalized_max_abs_error']<.12 # legacy blind spot
    assert metrics['OH_Y_trace_relative_max_error']>1000
    assert metrics['OH_Y_max_abs_difference']==pytest.approx(4.999e-9)
    assert metrics['OH_Y_active_points']==7


def test_no_reference_or_candidate_radicals_does_not_divide_by_zero():
    t=np.linspace(0,1,7)
    obj={'T_K':np.full(7,300.),'OH_Y':np.zeros(7),'H2O_Y':np.zeros(7),
         'ignition_observed':False,'threshold_crossing_time_s':None,
         'peak_dTdt_time_s':None,'max_time_step_s':1/6,
         'temperature_threshold_rise_K':20.,'ignition_time_s':None}
    assert _curve_error(obj,obj)['OH_Y_trace_relative_max_error']==0


def test_pinned_default_parameters_and_modified_case_rejection():
    plan=load_pinned_protocol('v10_h2o2_default')
    verify_protocol_parameters(plan,'h2o2.yaml',plan['training'],plan['holdout'],
        .006,80,(.4,.2,.1,.05),.1,.12)
    with pytest.raises(ValueError,match='differ'):
        verify_protocol_parameters(plan,'h2o2.yaml',plan['training'],plan['holdout'],
            .006,80,(.4,.2,.1,.05),.1,.9)
    with pytest.raises(ValueError,match='differ'):
        verify_protocol_parameters(plan,'h2o2.yaml',plan['training'],plan['holdout'],
            .1,80,(.4,.2,.1,.05),.1,.12)


def test_pinned_sweep_and_forged_sha_rejected():
    plan=load_pinned_protocol('v10_h2o2_composition')
    verify_protocol_parameters(plan,'h2o2.yaml',plan['training'],plan['holdout'],
        .006,80,(.4,.2,.1,.05),.1,.12)
    with pytest.raises(ValueError):
        load_pinned_protocol('anything_will_pass')
    with pytest.raises(ValueError):
        verify_summary_protocol(plan,{'protocol_id':plan['id'],
                                      'protocol_sha256':'f'*64})


def _fake_v10_summary(tmp_path,monkeypatch):
    result=sample_summary(tmp_path)
    result['evidence_schema']='v10-pinned'
    result['training_tolerance']=.12
    result['training_metrics']=[]
    for case in result['training_conditions']:
        row=dict(result['holdout_metrics'][0])
        row['condition_T_K_P_atm']=case
        row['OH_Y_trace_relative_max_error']=.02
        result['training_metrics'].append(row)
    result['training_max_normalized_error']=.03
    result['holdout_metrics'][0]['OH_Y_trace_relative_max_error']=.02
    # Make the controlled replay fixture represent a pinned, verified protocol.
    monkeypatch.setattr('plasmarom.protocol.load_pinned_protocol',
                        lambda *args: {'id':'fake','sha256':'stubbed'})
    monkeypatch.setattr('plasmarom.protocol.verify_summary_protocol',
                        lambda *args: None)
    monkeypatch.setattr(acceptance,'_load_mechanism_metadata',fake_verified_metadata)
    monkeypatch.setattr(acceptance,'_replay_training_observations',
                        lambda *args:[{k:v for k,v in row.items() if k!='condition_T_K_P_atm'}
                                      for row in result['training_ignition_observations']])
    monkeypatch.setattr(acceptance,'_replay_holdout_metrics',
                        lambda *args:[{k:v for k,v in row.items() if k!='condition_T_K_P_atm'}
                                      for row in (result['training_metrics'] if len(args[2])==2 else result['holdout_metrics'])])
    return result


def test_v10_training_candidate_replayed_independently(tmp_path,monkeypatch):
    report=_fake_v10_summary(tmp_path,monkeypatch)
    assert acceptance.assess_skeletal_summary(report,tmp_path)['accepted']
    replay=[{k:v for k,v in row.items() if k!='condition_T_K_P_atm'}
            for row in report['training_metrics']]
    replay[0]['OH_Y_trace_relative_max_error']=.8
    original=acceptance._replay_holdout_metrics
    monkeypatch.setattr(acceptance,'_replay_holdout_metrics',
                        lambda *args: replay if len(args[2])==2 else original(*args))
    decision=acceptance.assess_skeletal_summary(report,tmp_path)
    assert not decision['accepted']
    assert 'independent_training_candidate_replay_mismatch' in decision['reasons']


def test_v10_missing_training_candidate_evidence_fails_closed(tmp_path,monkeypatch):
    report=_fake_v10_summary(tmp_path,monkeypatch)
    del report['training_metrics']
    decision=acceptance.assess_skeletal_summary(report,tmp_path)
    assert not decision['accepted']
    assert 'missing_training_candidate_evidence' in decision['reasons']


def test_v10_unpinned_cannot_claim_acceptance(tmp_path,monkeypatch):
    report=_fake_v10_summary(tmp_path,monkeypatch)
    report['evidence_schema']='v10-unpinned-exploratory'
    decision=acceptance.assess_skeletal_summary(report,tmp_path)
    assert not decision['accepted']
    assert 'unlocked_exploratory_protocol_cannot_be_certified' in decision['reasons']

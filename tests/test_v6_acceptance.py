"""Independent regression tests of v6 fail-closed ignition and acceptance."""
import copy
import hashlib
import numpy as np
import pytest
from plasmarom.cantera_demo import _ignition_assessment, _curve_error, _max_curve_error, _safe_reaction_subset
from plasmarom.acceptance import assess_skeletal_summary


def test_no_temperature_rise_is_censored_not_ignition():
    t=np.linspace(0,.001,101)
    a=_ignition_assessment(t,300+0.1*t/t[-1])
    assert a['ignition_observed'] is False
    assert a['ignition_time_s'] is None
    assert a['ignition_status'].startswith('right_censored')


def test_intermediate_rise_is_ignition_proxy_not_true_ignition_claim():
    t=np.linspace(0,.01,101)
    T=300+200/(1+np.exp(-(t-.005)/.0005))
    result=_ignition_assessment(t,T)
    assert result['ignition_observed'] is True
    assert 0.003<result['ignition_time_s']<0.007


def test_endpoint_crossing_not_accepted():
    t=np.linspace(0,.01,101)
    T=300+21*(t/t[-1])**25
    r=_ignition_assessment(t,T)
    assert r['ignition_observed'] is False
    assert r['ignition_status']=='right_censored_crossing_near_end'


def _curve(ignition):
    t=np.linspace(0,.01,51)
    T=300+50/(1+np.exp(-(t-.005)/.0005)) if ignition else np.full_like(t,300.)
    d={'time_s':t,'T_K':T,'OH_Y':np.ones_like(t)*1e-7,
       'H2O_Y':np.ones_like(t)*1e-4}
    d.update(_ignition_assessment(t,T))
    return d


def test_nonignition_both_sides_not_validated():
    err=_curve_error(_curve(False),_curve(False))
    assert err['ignition_error_evaluable'] is False
    assert err['ignition_proxy_rel_error'] is None
    assert _max_curve_error(err)==0
    assert err['reference_ignition_observed'] is False


def test_different_ignition_status_is_penalized():
    err=_curve_error(_curve(True),_curve(False))
    assert err['ignition_proxy_rel_error']==1
    assert not err['ignition_classification_matches']


def valid_summary():
    # Reuse the stricter fixture; the old v6 fixture lacked provenance and
    # accepted a fake YAML just because its SHA-256 matched.
    # Deliberately represent the old report as incomplete to test rejection.
    return {'status':'HOLDOUT_PASS','full_species':10,'reduced_species':8,
            'full_reactions':30,'reduced_reactions':20,
            'accepted_with_ignition_evidence':True,
            'holdout_ignition_classifications_match':True,
            'ignition_training_cases_observed':1,
            'ignition_holdout_cases_observed':1,
            'holdout_max_normalized_error':.03,'holdout_tolerance':.12,
            'artifact':'holdout_pass_skeletal.yaml',
            'holdout_conditions':[[1500,1]],
            'holdout_metrics':[{'reference_ignition_observed':True,
                'candidate_ignition_observed':True,
                'ignition_classification_matches':True,
                'ignition_error_evaluable':True,
                'ignition_proxy_rel_error':.02,
                'T_K_normalized_max_abs_error':.03,
                'OH_Y_normalized_max_abs_error':.01,
                'H2O_Y_normalized_max_abs_error':.01}]}


def test_v6_incomplete_report_is_now_fail_closed(tmp_path):
    valid=valid_summary()
    blob=b'test-only fake artifact, not Cantera'
    (tmp_path/'holdout_pass_skeletal.yaml').write_bytes(blob)
    valid['artifact_sha256']=hashlib.sha256(blob).hexdigest()
    assessment=assess_skeletal_summary(valid,tmp_path)
    assert not assessment['accepted']
    assert 'invalid_training_conditions' in assessment['reasons']
    assert 'invalid_or_unreadable_cantera_mechanism' in assessment['reasons'] or 'cantera_unavailable_cannot_validate_mechanism' in assessment['reasons'] or 'missing_source_mechanism' in assessment['reasons']
    for key,value in [('status','NO REDUCTION'),('reduced_species',10),
                      ('holdout_max_normalized_error',float('nan')),
                      ('ignition_holdout_cases_observed',0),
                      ('holdout_ignition_classifications_match',False)]:
        bad=copy.deepcopy(valid);bad[key]=value
        assert not assess_skeletal_summary(bad,tmp_path)['accepted'],key


def test_censored_trajectories_do_not_pass_ignition_gate(tmp_path):
    obj=valid_summary();obj['ignition_holdout_cases_observed']=0
    obj['holdout_metrics'][0].update(reference_ignition_observed=False,
                                     candidate_ignition_observed=False,
                                     ignition_error_evaluable=False,
                                     ignition_proxy_rel_error=None)
    content=b'faux'
    (tmp_path/'holdout_pass_skeletal.yaml').write_bytes(content)
    obj['artifact_sha256']=hashlib.sha256(content).hexdigest()
    assert not assess_skeletal_summary(obj,tmp_path)['accepted']


class MockReaction:
    def __init__(self, reactants, products, efficiencies=None, collider='M'):
        self.reactants=reactants;self.products=products;self.orders={}
        if efficiencies is not None:
            self.third_body=type('T',(),{'collider':collider,'efficiencies':efficiencies})()
        else:self.third_body=None

class MockGas:
    def reactions(self):
        return [MockReaction({'H2':1},{'H':2}),
                MockReaction({'H':1},{'OH':1},{'AR':2}),
                MockReaction({'H':1},{'OH':1}, {'N2':.5})]


def test_third_body_reaction_not_silently_changed_by_pruning():
    rxns=_safe_reaction_subset(MockGas(),['H2','H','OH','N2'])
    assert len(rxns)==2
    assert 'AR' not in rxns[-1].third_body.efficiencies

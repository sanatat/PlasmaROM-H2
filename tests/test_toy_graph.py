import numpy as np
from plasmarom.kinetics import simulate
from plasmarom.toy_graph import NU,PARTICIPATION,graph_diagnostics


def test_stoichiometric_conservation():
    assert np.max(np.abs(NU.sum(axis=0)))==0
    assert PARTICIPATION.shape==(5,8)


def test_energy_budget_normalized_and_source_states_kept():
    result=graph_diagnostics(simulate(),threshold=.9,energy_tolerance=.01)
    assert abs(sum(result['toy_energy_fractions'].values())-1)<1e-10
    assert 'H2(v2) proxy' in result['energy_aware_candidate']

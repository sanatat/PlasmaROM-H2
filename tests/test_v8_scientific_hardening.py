"""Numerical and evidence-integrity regressions for version 8."""
from dataclasses import replace
from types import SimpleNamespace
import numpy as np
import pytest

from plasmarom.plasma_states import maxwell_reversible_rates
from plasmarom.electron_energy import ElectronEnergyConfig,run_coupled_electron_energy
from plasmarom.cantera_demo import _ignition_assessment,_curve_error
from plasmarom.acceptance import _reaction_fingerprint


def test_detailed_balance_coarse_energy_grids_and_thermal_conditions():
    gap=.55
    for spacing in (.7,.15):
        E=np.linspace(0,70,int(70/spacing)+1)
        sigma=np.where(E>gap,1e-20*(1-np.exp(-(E-gap)/.5)),0.)
        table=np.column_stack((E,sigma))
        for Te in (1.,2.,5.):
            up,down=maxwell_reversible_rates(table,Te,gap)
            assert up>0 and down>0
            assert down/up==pytest.approx(np.exp(gap/Te),rel=3e-15)


def test_degeneracies_and_bath_energy_remain_consistent():
    E=np.linspace(0,100,501)
    tbl=np.column_stack((E,np.where(E>.5,1e-21,0.)))
    a,b=maxwell_reversible_rates(tbl,2.,.5,degeneracy_lower=1,degeneracy_upper=3)
    assert b/a==pytest.approx(np.exp(.25)/3,rel=1e-13)
    with pytest.raises(ValueError,match='degeneracies'):
        maxwell_reversible_rates(tbl,2.,.5,degeneracy_upper=0)
    run=run_coupled_electron_energy(config=replace(ElectronEnergyConfig(),n_outputs=31))
    assert run['relative_energy_error']<1e-6
    assert run['min_population']>-1e-8


def test_threshold_crossing_and_slope_peak_are_disclosed_separately():
    time=np.arange(11,dtype=float)
    temperature=np.array([300,301,308,314,321,323,325,326,327,328,329.],float)
    result=_ignition_assessment(time,temperature)
    assert result['ignition_observed'] is True
    assert result['peak_dTdt_time_s']==pytest.approx(2.)
    assert result['threshold_crossing_time_s']==pytest.approx(3+6/7)
    assert result['ignition_time_s']==result['peak_dTdt_time_s']
    assert result['peak_dTdt_time_s']<result['threshold_crossing_time_s']
    assert result['max_time_step_s']==1.
    never=_ignition_assessment(time,np.full_like(time,300.))
    assert never['ignition_time_s'] is None
    assert never['threshold_crossing_time_s'] is None
    assert never['peak_dTdt_time_s'] is None


def test_ignition_clock_fields_are_preserved_in_curve_comparison():
    t=np.linspace(0,10,11)
    base=np.array([300,301,308,314,321,323,325,326,327,328,329.])
    def trajectory(T):
        d={'time_s':t,'T_K':T,'OH_Y':np.ones_like(t)*1e-7,
           'H2O_Y':np.ones_like(t)*1e-4}
        d.update(_ignition_assessment(t,T))
        return d
    metrics=_curve_error(trajectory(base),trajectory(base+0.5))
    assert metrics['reference_peak_dTdt_time_s']==2.
    assert metrics['reference_threshold_crossing_time_s']==pytest.approx(3+6/7)
    assert metrics['ignition_error_evaluable'] is True
    assert metrics['ignition_proxy_rel_error']==0


def test_full_rate_identity_changes_when_rate_or_third_body_changes():
    def make(rate=3.,efficiency=2.):
        rate_input={'equation':'H2 + M <=> H + H + M',
                    'rate-constant':{'A':rate,'b':0.,'Ea':0.},
                    'efficiencies':{'H2O':efficiency}}
        third=SimpleNamespace(collider='M',efficiencies={'H2O':efficiency},
                              default_efficiency=1.)
        return SimpleNamespace(equation='H2 + M <=> 2 H + M',
                               reactants={'H2':1.},products={'H':2.},
                               reversible=True,duplicate=False,orders={},
                               reaction_type='three-body',input_data=rate_input,
                               third_body=third)
    nominal=_reaction_fingerprint(make())
    assert nominal==_reaction_fingerprint(make())
    assert nominal!=_reaction_fingerprint(make(rate=3.1))
    assert nominal!=_reaction_fingerprint(make(efficiency=2.1))

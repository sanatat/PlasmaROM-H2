import numpy as np
import pytest
from plasmarom.verification import periodic_scalar_case, convergence_study
from plasmarom.plasma_states import (levels_eV,synthetic_table,maxwell_rate,
    integrate_levels,PlasmaConfig,read_cross_section_csv,run_plasma)
from plasmarom.thermal_flow import run_thermal
from plasmarom.energy_reduction import screen_power_channels
from plasmarom.hybrid_closure import lag_fraction,simulate_lag
from plasmarom.kinetics import Parameters
from plasmarom.cantera_demo import _safe_reaction_subset,run_skeletal


def test_convergence_expected_spatial_orders():
    data=convergence_study(ns=(16,32,64))
    assert .84<data['advection']['observed_order']<1.12
    assert 1.78<data['diffusion']['observed_order']<2.15


@pytest.mark.parametrize('mode',['advection','diffusion'])
def test_periodic_mass_conservation(mode):
    result=periodic_scalar_case(32,mode)
    assert result['mass_relative_error']<1e-12
    assert result['l2_error']>0


def test_real_h2_vibrational_energy_monotone():
    levels=levels_eV()
    assert levels[0]==0
    assert .4<levels[1]<.6
    assert levels[2]>levels[1]
    assert (levels[2]-levels[1])<(levels[1]-levels[0])


def test_maxwell_rate_constant_cross_section_limit():
    from scipy.constants import e,m_e
    te=2.
    energy=np.linspace(0,100,20001)
    sigma=np.full_like(energy,2e-20)
    numeric=maxwell_rate(np.column_stack([energy,sigma]),te,0.)
    expected=2e-20*np.sqrt(8*e*te/(np.pi*m_e))
    assert numeric==pytest.approx(expected,rel=2e-4)


def test_threshold_cutoff_makes_rate_smaller():
    E=np.linspace(0,60,1201);sigma=np.full_like(E,1e-20)
    rates=[maxwell_rate(np.column_stack((E,sigma)),2.,threshold) for threshold in (0.,.5,1.0,3.)]
    assert rates==sorted(rates,reverse=True)


def test_plasma_populations_and_energy_balance():
    levels=levels_eV();gap=np.diff(levels)
    sol=integrate_levels(synthetic_table(gap[0]),synthetic_table(gap[1]))
    assert sol['max_population_sum_error']<1e-8
    assert sol['relative_energy_closure']<1e-8
    assert sol['min_population']>-1e-9
    assert sol['fractions'][0,-1]<1
    assert sol['gas_heat_Jm3'][-1]>0


def test_plasma_zero_cross_section_remains_ground_state():
    table=np.column_stack([np.linspace(0,50,21),np.zeros(21)])
    sol=integrate_levels(table,table)
    assert np.allclose(sol['fractions'][0],1)
    assert np.allclose(sol['fractions'][1:],0)
    assert np.allclose(sol['electron_energy_Jm3'],0)


def test_invalid_cross_sections_and_conditions():
    tab=synthetic_table(.5)
    with pytest.raises(ValueError):maxwell_rate(tab,-1,.5)
    with pytest.raises(ValueError):maxwell_rate(tab,2,-.5)
    tab[6,1]=-1
    with pytest.raises(ValueError):maxwell_rate(tab,2,.5)
    with pytest.raises(ValueError):integrate_levels(synthetic_table(.5),synthetic_table(.5),PlasmaConfig(ne_m3=-1))


def test_cross_section_csv_provenance_required(tmp_path):
    csv=tmp_path/'cross.csv'
    csv.write_text('energy_eV,sigma_m2\n0,0\n1,0\n2,1e-20\n')
    arr=read_cross_section_csv(csv)
    assert arr.shape==(3,2)
    with pytest.raises(ValueError):run_plasma(tmp_path/'out',csv,csv,source01='dummy source')
    bad=tmp_path/'invalid.csv'
    bad.write_text('energy_eV,sigma_m2\n0,0\n2,0\n1,0\n')
    with pytest.raises(ValueError):read_cross_section_csv(bad)


def test_vibration_demo_is_explicitly_synthetic(tmp_path):
    result=run_plasma(tmp_path)
    assert 'SYNTHETIC' in result['data_status']
    assert (tmp_path/'plasma_results.json').exists()
    assert (tmp_path/'plasma_states.png').exists()


def test_conservative_dimensionless_heat(tmp_path):
    row=run_thermal(tmp_path,n=16,dt=.015,duration=.3)
    assert row['max_global_energy_proxy_drift']<1e-10
    assert row['max_total_species_drift']<1e-10
    assert row['min_species']>=-1e-12
    assert (tmp_path/'thermal_surrogate.png').exists()


def test_invalid_thermal_setup():
    with pytest.raises(ValueError):run_thermal(n=17)
    with pytest.raises(ValueError):run_thermal(heat_release=-1)


def test_energy_screen_preserves_target_share():
    t=np.linspace(0,1,30)
    rates=np.repeat(np.array([[10.,2.,.1]]),len(t),axis=0)
    result=screen_power_channels(t,rates,[1,3,1],['A','B','C'],fraction=.9)
    assert result['retained_fraction']>=.9
    assert result['selected']==['A','B']


def test_energy_screen_invalid_inputs():
    with pytest.raises(ValueError):screen_power_channels([0,1],[[1],[2]],[1],['A'],fraction=1.2)
    with pytest.raises(ValueError):screen_power_channels([0,1],[[1],[-2]],[1],['A'])


def test_lag_model_conservation():
    p=Parameters()
    assert 0<=lag_fraction(1.2,p,.1)<=1
    sol=simulate_lag(p,tau=.02)
    assert sol.success
    assert np.max(abs(sol.y.sum(axis=0)-1))<1e-8
    assert sol.y.min()>-1e-8
    with pytest.raises(ValueError):simulate_lag(p,tau=-.01)


def test_cantera_skeletal_parameter_validation_without_cantera(tmp_path):
    for kw in ({'training':()},{'holdout':()},{'n_snapshots':1},
               {'thresholds':(2.,)},{'t_end':0},{'train_tolerance':-1}):
        with pytest.raises(ValueError):run_skeletal(tmp_path,**kw)


def test_third_body_dependency_filter_no_cantera():
    class Third:
        collider='M'
        efficiencies={'N2':2.,'AR':1.5}
    class Reaction:
        reactants={'H':1.,'O2':1.}
        products={'HO2':1.}
        orders={}
        third_body=Third()
    class Phase:
        def reactions(self):return [Reaction()]
    assert len(_safe_reaction_subset(Phase(),['H','O2','HO2','N2']))==0
    assert len(_safe_reaction_subset(Phase(),['H','O2','HO2','N2','AR']))==1

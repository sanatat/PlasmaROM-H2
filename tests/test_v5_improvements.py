"""Regression + falsification tests added for version 5."""
from dataclasses import replace
import numpy as np
import pytest
from scipy.integrate import trapezoid

from plasmarom.plasma_states import maxwell_rate, synthetic_table, levels_eV
from plasmarom.eedf import validate_eedf, rate_from_eedf, maxwellian_pdf
from plasmarom.electron_energy import ElectronEnergyConfig, run_coupled_electron_energy, run_energy_demo
from plasmarom.channel_reduction import evaluate_channel_pruning


def test_maxwell_tail_coverage_rejects_short_tables():
    table=synthetic_table(.5)
    short=table[table[:,0]<=5]
    with pytest.raises(ValueError,match='Insufficient EEDF'):
        maxwell_rate(short,te_eV=2.,threshold_eV=.5)
    assert maxwell_rate(table,2.,.5)>0


def test_maxwellian_eedf_convolution_matches_original_rate():
    E=np.linspace(0,80,16001)
    pdf=np.column_stack((E,maxwellian_pdf(E,2.)))
    sig=synthetic_table(.5)
    assert validate_eedf(pdf)==pytest.approx(1.,rel=1e-3)
    assert rate_from_eedf(sig,pdf,.5)==pytest.approx(maxwell_rate(sig,2.,.5),rel=2e-4)


def test_unnormalized_pdf_fails_and_unknown_high_energy_fails():
    E=np.linspace(0,80,801)
    p=np.column_stack((E,maxwellian_pdf(E,2.)))
    with pytest.raises(ValueError,match='unity'):
        validate_eedf(np.column_stack((E,p[:,1]*1.5)))
    E2=np.linspace(0,300,6001)
    extended=np.column_stack((E2,maxwellian_pdf(E2,20.)))
    with pytest.raises(ValueError,match='outside cross-section'):
        rate_from_eedf(synthetic_table(.5),extended,.5)


@pytest.mark.parametrize('power',[0.,100.,3000.,6000.])
def test_coupled_electron_energy_and_populations_are_conserved(power):
    result=run_coupled_electron_energy(config=ElectronEnergyConfig(input_power_W_m3=power))
    assert result['relative_energy_error']<1e-6
    assert result['max_population_sum_error']<1e-8
    assert result['min_population']>-1e-8
    assert np.all(np.isfinite(result['te_eV']))
    assert np.all(result['te_eV']>0)
    assert result['injected_energy_J_m3'][-1]==pytest.approx(power*2e-4,rel=1e-8,abs=1e-8)


def test_plasma_is_not_a_fixed_temperature_accounting_identity():
    cfg=ElectronEnergyConfig(input_power_W_m3=0.)
    d=run_coupled_electron_energy(config=cfg)
    assert d['te_eV'][-1] < cfg.te_initial_eV*0.1
    excited=run_coupled_electron_energy(config=replace(cfg,input_power_W_m3=6000.))
    assert excited['te_eV'][-1]>d['te_eV'][-1]


def test_external_eedf_and_cross_section_inputs_must_be_paired():
    with pytest.raises(ValueError):
        run_coupled_electron_energy(synthetic_table(.5),None)
    with pytest.raises(ValueError):
        run_coupled_electron_energy(config=ElectronEnergyConfig(ne_m3=-1))


def test_synthetic_channel_pruning_meets_heldout_error(tmp_path):
    d=evaluate_channel_pruning(tmp_path)
    assert d['status']=='HELDOUT_PASS'
    assert d['removed_transitions']==['v1<->v2']
    assert d['retained_energy_fraction']>=d['energy_selection_fraction_target']
    assert d['worst_holdout_error']<=d['allowed_max_normalized_error']
    assert (tmp_path/'channel_reduction_summary.json').exists()


def test_requested_unrealistic_energy_fraction_refuses_pruning(tmp_path):
    d=evaluate_channel_pruning(tmp_path, energy_fraction=1.)
    assert d['status']=='NOT_VALIDATED_REDUCTION'
    assert d['removed_transitions']==[]


def test_plasma_energy_report_round_trip(tmp_path):
    d=run_energy_demo(tmp_path)
    assert d['relative_energy_error']<1e-6
    assert (tmp_path/'electron_energy_summary.json').exists()
    assert (tmp_path/'electron_energy.png').exists()


def test_gross_flux_screen_protects_reversible_exchange():
    from plasmarom.reduction import (direct_interaction_coefficients,
                                      gross_flux_interaction_coefficients)
    nu=np.array([[1.],[-1.],[0.]])
    participation=np.array([[True],[True],[True]])
    net=direct_interaction_coefficients(nu,np.array([0.]),participation)
    gross=gross_flux_interaction_coefficients(nu,np.array([10.]),np.array([10.]),participation)
    assert net[0,1]==0.
    assert gross[0,1]==pytest.approx(1.)
    with pytest.raises(ValueError):
        gross_flux_interaction_coefficients(nu,np.array([-1.]),np.array([2.]),participation)


def test_cantera_flame_input_checks_without_cantera(tmp_path):
    from plasmarom.cantera_flame import run_flame
    with pytest.raises(ValueError):run_flame(tmp_path,width_m=-1)
    with pytest.raises(ValueError):run_flame(tmp_path,initial_T_K=0)


def test_user_eedf_rate_requires_source(tmp_path):
    from plasmarom.eedf import run_eedf_rate
    with pytest.raises(ValueError,match='sources'):
        run_eedf_rate(tmp_path,'missing.csv','missing.csv','','',.5)


def test_local_multicore_results_match_serial(tmp_path):
    from plasmarom.parallel import compare_parallel
    row=compare_parallel(tmp_path,powers=(500.,3000.),workers=2)
    assert row['max_absolute_serial_parallel_difference']==0
    assert row['serial_wall_s']>0 and row['parallel_wall_s']>0


def test_sample_external_eedf_csv_and_provenance(tmp_path):
    from pathlib import Path
    from plasmarom.eedf import run_eedf_rate
    base=Path(__file__).resolve().parents[1]/'examples'/'data'
    result=run_eedf_rate(tmp_path,
        base/'SYNTHETIC_Maxwellian_EEDF_Te2eV.csv',
        base/'SYNTHETIC_H2_v0_v1_cross_section.csv',
        'SYNTHETIC EEDF', 'SYNTHETIC cross section', levels_eV()[1])
    assert result['rate_coefficient_m3_s']>0
    assert len(result['eedf_sha256'])==64
    assert (tmp_path/'external_eedf_rate.json').exists()

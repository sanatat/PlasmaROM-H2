import csv
import json
import numpy as np
import pytest
from plasmarom.ammonia_h2 import check_initial_mixture
from plasmarom.energy_graph import energy_aware_screen
from plasmarom.experiment_compare import evaluate_experiment
from plasmarom.reactive_transport_1d import simulate
from plasmarom.bolos_bridge import run_bolos

def test_nh3_h2_stoichiometric_phi():
    assert check_initial_mixture({'NH3':1,'H2':1,'O2':1.25,'N2':4.7})==pytest.approx(1)

def test_nh3_h2_requires_oxidizer():
    with pytest.raises(ValueError):check_initial_mixture({'NH3':1,'H2':1,'O2':0,'N2':5})

def test_energy_graph_keeps_power_dominant_channel():
    g=np.array([[0.,.1,0],[0,0,.05],[0,0,0]])
    r=np.tile([1.,50.,1.],(10,1))
    result=energy_aware_screen(g,['v0'],r,[1,5,3],np.linspace(0,1,10),['v0','v1','v2'],graph_threshold=.4,min_energy_fraction=.98)
    assert 'v0' in result['selected'] and 'v1' in result['selected']
    assert result['retained_energy_fraction']>=.98
    assert result['status'].startswith('SCREENING_ONLY')

def test_energy_graph_rejects_nan():
    with pytest.raises(ValueError):energy_aware_screen(np.full((2,2),np.nan),['A'],np.ones((3,2)),[1,2],[0,.5,1],['A','B'])

def test_1d_transport_conserves_mass_and_energy():
    d=simulate(cells=80,end_s=.012)
    assert d['summary']['global_mass_drift']<2e-13
    assert d['summary']['global_total_enthalpy_drift_J_kg']<1e-7
    assert np.all(d['A']>=0) and np.all(d['B']>=0)

def test_uniform_reaction_against_analytic():
    end=.025;rate=7.;a0=.6
    d=simulate(cells=24,end_s=end,k_ref_s=rate,activation_J_mol=0,uniform=True,u_m_s=0,diffusivity_m2_s=0,thermal_diffusivity_m2_s=0)
    expected=a0*np.exp(-rate*end)
    assert max(abs(d['A']-expected))<1e-10
    assert max(abs(d['T_K']-(600+2e5/1000*(a0-expected))))<1e-8

def test_experiment_audit_requires_real_provenance(tmp_path):
    a=tmp_path/'measured.csv';b=tmp_path/'predicted.csv'
    with a.open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['case_id','observable','measured_value','uncertainty_1sigma','units','source_doi_or_url','data_origin']);
        w.writerow(['case1','ignition_s',.02,.001,'s','','SYNTHETIC'])
    with b.open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['case_id','observable','predicted_value','units']);w.writerow(['case1','ignition_s',.022,'s'])
    with pytest.raises(ValueError):evaluate_experiment(a,b,tmp_path/'result.json')

def test_experiment_audit_rejects_mismatched_units(tmp_path):
    a=tmp_path/'m.csv';b=tmp_path/'p.csv'
    a.write_text('case_id,observable,measured_value,uncertainty_1sigma,units,source_doi_or_url,data_origin\nc,flame_speed,2,0.1,m/s,https://doi.org/example,MEASURED_INDEPENDENT\n')
    b.write_text('case_id,observable,predicted_value,units\nc,flame_speed,2,cm/s\n')
    with pytest.raises(ValueError):evaluate_experiment(a,b,tmp_path/'a.json')

def test_bolos_bridge_fails_on_missing_collisions(tmp_path):
    with pytest.raises(FileNotFoundError):run_bolos(tmp_path/'no.xs',tmp_path/'output',source_url='https://example.org/collisions',mixture={'H2':1.0})

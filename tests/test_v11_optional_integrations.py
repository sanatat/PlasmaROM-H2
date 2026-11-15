"""Real numerical integration tests; skip only when Cantera is unavailable.
CI must import Cantera FIRST, therefore it cannot silently skip these tests.
"""
import pytest
ct=pytest.importorskip('cantera')

def test_real_nh3_h2_mechanism_available():
    gas=ct.Solution('example_data/ammonia-CO-H2-Alzueta-2023.yaml','baseline')
    assert {'NH3','H2','O2','N2','NO','OH'}<=set(gas.species_names)
    assert gas.n_reactions>100

def test_real_boltzmann_eedf_from_phelps(tmp_path):
    from plasmarom.cantera_boltzmann import run_phelps_eedf
    data=run_phelps_eedf(tmp_path,fields_Td=(200.,))
    assert data['results'][0]['mean_energy_eV']>0
    assert data['status'].startswith('REAL_CANTERA_TWO_TERM')

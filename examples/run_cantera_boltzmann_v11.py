"""Cantera 3.2 real two-term EEDF on published Phelps N2/O2 example."""
from plasmarom.cantera_boltzmann import run_phelps_eedf
if __name__=='__main__':
    print(run_phelps_eedf('results_v11/cantera_boltzmann'))

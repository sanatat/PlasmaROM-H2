"""Runs Cantera 3.2 Alzueta et al. NH3/H2 reactor; no fabricated outputs."""
from plasmarom.ammonia_h2 import run_ammonia_h2
if __name__=='__main__':
    result=run_ammonia_h2('results_v11/ammonia_h2')
    print(result['status'], result['cases'])

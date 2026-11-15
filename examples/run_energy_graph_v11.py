"""Reproducible synthetic energy-aware graph screen; NOT published P-DRGEP."""
from pathlib import Path
import json
import numpy as np
from plasmarom.energy_graph import energy_aware_screen
if __name__=='__main__':
    labels=['H2_v0','H2_v1','H2_v2','dissociation']
    graph=np.array([[0.,.7,0.,.05],[.2,0.,.2,.1],[0.,.1,0.,.5],[0,0,.1,0.]])
    t=np.linspace(0,1e-4,101)
    rates=np.tile([30e18,10e18,1e18,1e16],(t.size,1))*np.exp(-t[:,None]/1e-4)
    energy=[.5,.5,.5,4.48]  # Demonstration values; NOT cross-section-derived
    result=energy_aware_screen(graph,['H2_v0'],rates,energy,t,labels,min_energy_fraction=.95)
    out=Path('results_v11');out.mkdir(exist_ok=True)
    (out/'SYNTHETIC_energy_graph.json').write_text(json.dumps(result,indent=2))
    print(result['status'],result['selected'],result['retained_energy_fraction'])

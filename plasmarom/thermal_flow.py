"""Dimensionless 2D reaction-transport + heat-release proxy.

A conservative scalar Tproxy is raised by H * local increment of product P.
The combined potential + sensible energy proxy Tproxy+H*(1-P) is globally
constant with periodic boundaries and matched FV transport. This is NOT
heat capacity, enthalpy chemistry, variable-density combustion, or a flame.
"""
from pathlib import Path
import json
import numpy as np
from .flow import (FlowConfig,SpectralFlow,divergence_free_face_velocities,
                   finite_volume_cfl,scalar_transport,reaction_ssprk2)
from .kinetics import Parameters


def run_thermal(out=None, n=32, dt=0.014, duration=1.4, heat_release=1.0,
                diffusion=0.005, parameters=Parameters()):
    if n<12 or n%2 or dt<=0 or duration<=0 or heat_release<=0 or diffusion<0:
        raise ValueError('Invalid thermal demonstration settings')
    flow=SpectralFlow(n,nu=0.01)
    mask=np.exp(-((flow.x-np.pi)**2+(flow.y-np.pi)**2)/(2*0.9**2))
    y=np.zeros((5,n,n));y[0]=1
    # Initial sensible scalar is zero; chemical potential is H*(1-P).
    theta=np.zeros((n,n))
    initial_energy=float(np.mean(theta+heat_release*(1-y[-1])))
    t=0; peak_balance=0;min_species=1
    for j in range(int(np.ceil(duration/dt))):
        hstep=min(dt,duration-t)
        if hstep<=0:break
        faces=divergence_free_face_velocities(flow.u,flow.v,flow.dx)
        if hstep*finite_volume_cfl(faces,flow.dx,diffusion)>.85:
            raise ValueError('Transport CFL exceeded')
        p0=y[-1].copy()
        y,_=reaction_ssprk2(y,t,hstep/2,parameters,mask)
        theta+=heat_release*(y[-1]-p0)
        y+=hstep*scalar_transport(y,flow.u,flow.v,flow.dx,diffusion,faces)
        theta+=hstep*scalar_transport(theta,flow.u,flow.v,flow.dx,diffusion,faces)
        p0=y[-1].copy()
        y,_=reaction_ssprk2(y,t+hstep/2,hstep/2,parameters,mask)
        theta+=heat_release*(y[-1]-p0)
        flow.step(hstep)
        t+=hstep
        peak_balance=max(peak_balance,abs(float(np.mean(theta+heat_release*(1-y[-1])))-initial_energy))
        min_species=min(min_species,float(y.min()))
    result={'model':'dimensionless coupled heat proxy, NOT thermochemical combustion',
            'n':n,'t_end':t,'min_species':min_species,'min_temperature_proxy':float(theta.min()),
            'max_temperature_proxy':float(theta.max()),'max_global_energy_proxy_drift':peak_balance,
            'max_total_species_drift':float(np.max(abs(y.sum(axis=0)-1)))}
    if out is not None:
        out=Path(out);out.mkdir(parents=True,exist_ok=True)
        (out/'thermal_results.json').write_text(json.dumps(result,indent=2))
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig,axes=plt.subplots(1,2,figsize=(8.4,3.6))
        for ax,z,title in zip(axes,(y[-1],theta),('Product-equivalent pool','Sensible-energy proxy')):
            im=ax.imshow(z.T,origin='lower',extent=(0,2*np.pi,0,2*np.pi))
            ax.set_title(title);fig.colorbar(im,ax=ax,shrink=.76)
        fig.suptitle('Illustrative energy-release transport — NOT flame CFD')
        fig.tight_layout();fig.savefig(out/'thermal_surrogate.png',dpi=150);plt.close(fig)
    return result

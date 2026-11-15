"""Manufactured smooth periodic scalar advection/diffusion checks.

Tests distinguish consistency with a known PDE solution from physical
validation of combustion. Velocity is held fixed and divergence free.
"""
from pathlib import Path
import json
import numpy as np
from .flow import scalar_transport, finite_volume_cfl, divergence_free_face_velocities


def periodic_scalar_case(n=32, mode='advection', t_end=0.35):
    if n<8 or n%2 or mode not in ('advection','diffusion') or t_end<=0:
        raise ValueError('Invalid grid/mode/time')
    dx=2*np.pi/n
    x=dx*np.arange(n)
    xx,yy=np.meshgrid(x,x,indexing='ij')
    speed=1.0 if mode=='advection' else 0.0
    diffusion=0.12 if mode=='diffusion' else 0.0
    u=np.full((n,n),speed);v=np.zeros_like(u)
    faces=divergence_free_face_velocities(u,v,dx)
    max_dt=0.35/max(finite_volume_cfl(faces,dx,diffusion),1e-12)
    # Reduce temporal truncation sufficiently to expose the spatial scheme's order.
    dt=min(max_dt,0.015*dx if mode=='advection' else 0.02*dx**3)
    nsteps=int(np.ceil(t_end/dt))
    dt=t_end/nsteps
    state=0.5+0.22*np.sin(xx)*np.cos(yy)
    mass_initial=float(state.mean())
    for _ in range(nsteps):
        state+=dt*scalar_transport(state,u,v,dx,diffusion,faces)
    exact=0.5+0.22*np.sin(xx-speed*t_end)*np.cos(yy)*np.exp(-2*diffusion*t_end)
    return {'n':n,'mode':mode,'t_end':t_end,'steps':nsteps,'dt':dt,
            'l2_error':float(np.sqrt(np.mean((state-exact)**2))),
            'max_abs_error':float(np.max(np.abs(state-exact))),
            'mass_relative_error':float(abs(state.mean()-mass_initial)/mass_initial)}


def convergence_study(ns=(16,32,64), out=None):
    data={}
    for mode in ('advection','diffusion'):
        rows=[periodic_scalar_case(n,mode) for n in ns]
        x=np.log(np.asarray([r['n'] for r in rows],float))
        y=np.log(np.asarray([r['l2_error'] for r in rows]))
        order=float(-np.polyfit(x,y,1)[0])
        data[mode]={'rows':rows,'observed_order':order}
    if out is not None:
        p=Path(out);p.mkdir(parents=True,exist_ok=True)
        (p/'numerical_verification.json').write_text(json.dumps(data,indent=2))
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig,ax=plt.subplots(figsize=(6,4))
        for mode,c in data.items():
            ax.loglog([r['n'] for r in c['rows']], [r['l2_error'] for r in c['rows']], 'o-',
                      label=f'{mode}: order {c["observed_order"]:.2f}')
        ax.set(xlabel='Grid width N (N×N)',ylabel='L2 error',title='Periodic scalar spatial verification')
        ax.grid(True,which='both',alpha=.25);ax.legend();fig.tight_layout()
        fig.savefig(p/'spatial_convergence.png',dpi=160);plt.close(fig)
    return data

"""Quantitative verification and plots for an educational surrogate model."""
import json
import time
from pathlib import Path
import numpy as np
from .kinetics import Parameters, simulate, full_observables, pca_variance, FULL_NAMES, reaction_flux_screen
from .flow import FlowConfig, run_flow
from .toy_graph import graph_diagnostics
from .learned_closure import fit_closure, held_out_metrics


def zero_d(out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    p=Parameters()
    full=simulate(p)
    reduced=simulate(p,reduced=True)
    if not full.success or not reduced.success:raise RuntimeError('Integration failure')
    obs=full_observables(full)
    err=np.max(np.abs(obs-reduced.y),axis=1)
    from .kinetics import qssa_reconstruct
    reconstructed=qssa_reconstruct(reduced,p)
    state_errors=np.max(np.abs(full.y-reconstructed),axis=1)
    cases=[]
    for strength in (0.5,1.,1.5):
        for width in (0.2,0.34,0.55):
            q=Parameters(pulse_strength=strength,pulse_width=width)
            a,b=simulate(q),simulate(q,reduced=True)
            cases.append({'strength':strength,'width':width,
                          'MAE_R':float(np.mean(abs(a.y[3]-b.y[2]))),
                          'max_abs_R':float(np.max(abs(a.y[3]-b.y[2]))),
                          'MAE_P':float(np.mean(abs(a.y[4]-b.y[3]))),
                          'max_abs_P':float(np.max(abs(a.y[4]-b.y[3])))})
    times={}
    for kind in ('full','QSSA'):
        runs=[]
        for _ in range(8):
            start=time.perf_counter()
            simulate(p,reduced=(kind=='QSSA'),n_points=201)
            runs.append(time.perf_counter()-start)
        times[kind+'_seconds_median']=float(np.median(runs))
    ratio,axes=pca_variance(np.vstack([simulate(Parameters(pulse_strength=s)).y.T for s in (0.5,1,1.5)]))
    fig,ax=plt.subplots(figsize=(8,4.6))
    for i in (0,1,2,3,4):ax.plot(full.t,full.y[i],label=FULL_NAMES[i])
    ax.plot(reduced.t,reduced.y[3],'k--',lw=1.8,label='QSSA converted pool')
    ax.set(xlabel='Dimensionless time',ylabel='Hydrogen-equivalent fraction',title='State-resolved surrogate kinetics (NOT calibrated plasma chemistry)')
    ax.legend(fontsize=7,ncol=2);fig.tight_layout();fig.savefig(out/'kinetics.png',dpi=150);plt.close(fig)
    graph=graph_diagnostics(full,p)
    fit_weights,fit_residual=fit_closure()
    heldout=held_out_metrics(fit_weights)
    channels=reaction_flux_screen(full.t,full.y,p)
    fig,ax=plt.subplots(figsize=(8,4.3))
    ax.barh([n for n,v in reversed(channels)],[v for n,v in reversed(channels)])
    ax.set(xlabel='Integrated transition flux (model units)',title='Reaction-channel screening (not P-DRGEP)')
    fig.tight_layout();fig.savefig(out/'pathway_flux.png',dpi=150);plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4))
    ax.bar(np.arange(1,len(ratio)+1),np.cumsum(ratio)*100)
    ax.set(xlabel='Number of PCA components',ylabel='Cumulative explained variance (%)',ylim=(0,105),title='PCA diagnostic (not a kinetic solver)')
    fig.tight_layout();fig.savefig(out/'pca.png',dpi=150);plt.close(fig)
    summary={'model':'dimensionless illustrative hydrogen-equivalent excitation/oxidation; not real PAC',
             'full_species':5,'qssa_species':4,'species_reduction_fraction':0.2,
             'max_abs_error_F0_B_R_P':err.tolist(),
             'max_abs_reconstructed_state_errors_F0_F1_F2_R_P':state_errors.tolist(),
             'mass_residual_full':float(np.max(abs(full.y.sum(axis=0)-1))),
             'mass_residual_QSSA':float(np.max(abs(reduced.y.sum(axis=0)-1))),
             'min_population_full':float(full.y.min()),
             'pca_explained_variance':ratio.tolist(),
             'integrated_pathway_flux_ranking':channels,
             'toy_drgep_and_energy_target_screening':graph,
             'learned_nonnegative_closure':{'basis':'drive and drive squared',
                 'weights':fit_weights.tolist(),'training_residual_l2':fit_residual,
                 'training_cases':6,'heldout_cases':heldout,
                 'warning':'Toy NNLS regression, not physics-informed neural networks or calibrated plasma kinetics.'},
             'max_omitted_F2_population_full':float(full.y[2].max()),
             'validation_grid':cases,
             'runtimes':times,
             'runtime_warning':'Small systems may not speed up under Python/SciPy overhead; do not infer industrial speedup.'}
    (out/'kinetics_results.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    return summary


def spatial(out, quick=True):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    fine_n,coarse_n = (32,16) if quick else (64,32)
    fine_dt=0.025 if quick else 0.0125
    configurations = [
        (fine_n,fine_dt,0.,False),
        (fine_n,fine_dt,0.,True),
        (coarse_n,2*fine_dt,0.,True),
        (coarse_n,2*fine_dt,0.18,True),
    ]
    data=[]
    for n,dt,cs,reduced in configurations:
        result=run_flow(FlowConfig(n=n,dt=dt,t_end=4.2,cs=cs,reduced=reduced))
        data.append(result)
    titles=('Full / fine grid','QSSA / fine grid','QSSA / coarse grid','QSSA / coarse LES-style closure')
    fig,axs=plt.subplots(1,4,figsize=(16,3.8))
    high=max(d['species'][-1].max() for d in data)
    for ax,d,title in zip(axs,data,titles):
        z=d['species'][-1]
        im=ax.imshow(z.T,origin='lower',extent=(0,2*np.pi,0,2*np.pi),vmin=0,vmax=high)
        ax.set(xlabel='x',ylabel='y',title=title)
    fig.colorbar(im,ax=axs,shrink=0.78,label='Converted-equivalent fraction')
    fig.suptitle('2D Navier-Stokes / surrogate chemical transport (no real flame or validated LES)')
    fig.subplots_adjust(left=0.06,right=0.91,bottom=0.15,top=0.78,wspace=0.28)
    fig.savefig(out/'reactive_flow.png',dpi=150);plt.close(fig)
    summary=[]
    for d in data:
        summary.append({k:float(d[k]) if k not in ('n','reduced') else d[k]
                        for k in ('n','t','reduced','cs','peak_divergence','min_pool','sum_deviation','taylor_green_l2','max_chemical_substeps_per_half_step')})
    full_fine=data[0]['species'][-1]
    summary[1]['converted_pool_MAE_vs_full_fine']=float(np.mean(abs(data[1]['species'][-1]-full_fine)))
    for i in (2,3):
        baseline_downsampled=full_fine[::2,::2]
        summary[i]['converted_pool_MAE_vs_full_fine_downsampled']=float(np.mean(abs(data[i]['species'][-1]-baseline_downsampled)))
    # Same-n grid with halved time step: sensitivity diagnostic, not an
    # independent convergence-order proof (spatial error is unchanged).
    check=run_flow(FlowConfig(n=16,dt=.025,t_end=1.2))
    refined=run_flow(FlowConfig(n=16,dt=.0125,t_end=1.2))
    summary.append({'comparison':'same-grid temporal refinement',
                    'n':16,'coarse_dt':.025,'fine_dt':.0125,
                    'product_MAE':float(np.mean(abs(check['species'][-1]-refined['species'][-1])))})
    from .kinetics import Parameters
    strong=run_flow(FlowConfig(n=20,dt=.02,t_end=1.5),Parameters(pulse_strength=100.))
    summary.append({'comparison':'high-excitation numerical stress test (synthetic)',
                    'pulse_strength':100.,'n':20,'dt':.02,
                    'min_pool':strong['min_pool'],
                    'sum_deviation':strong['sum_deviation'],
                    'max_substeps':strong['max_chemical_substeps_per_half_step']})
    (out/'flow_results.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    return summary

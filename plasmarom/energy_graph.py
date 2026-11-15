"""Auditable energy-aware graph screening research prototype.

NOT the published plasma P-DRGEP algorithm. It demonstrates a constrained
union of graph importance and electron-energy channel retention on given data.
"""
from __future__ import annotations
import numpy as np
from scipy.integrate import trapezoid


def energy_aware_screen(interaction, targets, channel_rates, energy_losses_eV,
                        t_s, channel_names, *, graph_threshold=0.1,
                        min_energy_fraction=0.99):
    g=np.asarray(interaction,float)
    r=np.asarray(channel_rates,float)
    e=np.asarray(energy_losses_eV,float)
    t=np.asarray(t_s,float)
    n=len(channel_names)
    if (g.shape!=(n,n) or len(set(channel_names))!=n or r.shape!=(t.size,n) or
        e.shape!=(n,) or t.ndim!=1 or t.size<2 or
        np.any(~np.isfinite(g)) or np.any(g<0) or np.any(g>1) or
        np.any(~np.isfinite(r)) or np.any(r<0) or np.any(~np.isfinite(e)) or
        not np.isfinite(t).all() or np.any(np.diff(t)<=0) or
        not 0<min_energy_fraction<=1 or not 0<=graph_threshold<=1):
        raise ValueError('Invalid graph/rate/energy data')
    if not targets or any(x not in channel_names for x in targets):
        raise ValueError('Invalid target nodes')
    importance=np.zeros(n)
    # widest-path propagation; max-product graph over all target nodes
    for source in targets:
        weights=np.zeros(n);weights[channel_names.index(source)]=1
        for _ in range(n):
            next_weights=np.maximum(weights,np.max(weights[:,None]*g,axis=0))
            if np.allclose(next_weights,weights,atol=0,rtol=0): break
            weights=next_weights
        importance=np.maximum(importance,weights)
    # Signed energy exchange: positive => electron energy lost to the gas;
    # negative => electron energy gained through superelastic collisions.
    # Do not allow the two directions to cancel during screening.
    energy=trapezoid(r*e[None,:],t,axis=0)
    losses=np.maximum(energy,0.)
    gains=np.maximum(-energy,0.)
    loss_total=float(losses.sum())
    gain_total=float(gains.sum())
    selected={name for name,score in zip(channel_names,importance) if score>=graph_threshold}
    for magnitudes,total in ((losses,loss_total),(gains,gain_total)):
        if total<=0:
            continue
        retained=sum(magnitudes[channel_names.index(name)] for name in selected)
        for idx in np.argsort(-magnitudes,kind='stable'):
            if retained/total>=min_energy_fraction-1e-12:
                break
            if magnitudes[idx]<=0:
                continue
            selected.add(channel_names[idx])
            retained=sum(magnitudes[channel_names.index(name)] for name in selected)
    selected_indices=[i for i,name in enumerate(channel_names) if name in selected]
    retained_loss=float(losses[selected_indices].sum())
    retained_gain=float(gains[selected_indices].sum())
    absolute_total=loss_total+gain_total
    return {
        'selected':[name for name in channel_names if name in selected],
        'graph_importance':dict(zip(channel_names,importance.tolist())),
        'integrated_electron_energy_eV_units':dict(zip(channel_names,energy.tolist())),
        'integrated_energy_loss_eV_units':dict(zip(channel_names,losses.tolist())),
        'integrated_energy_gain_eV_units':dict(zip(channel_names,gains.tolist())),
        'retained_loss_fraction':float(retained_loss/loss_total) if loss_total>0 else None,
        'retained_gain_fraction':float(retained_gain/gain_total) if gain_total>0 else None,
        'retained_energy_fraction':float((retained_loss+retained_gain)/absolute_total) if absolute_total>0 else None,
        'net_energy_exchange_eV_units':float(energy.sum()),
        'energy_sign_convention':'positive=electron energy loss; negative=electron energy gain (superelastic)',
        'status':'SCREENING_ONLY_NOT_VALIDATED_P_DRGEP'
    }

"""Directed Relation Graph with Error Propagation (DRGEP) core.

A real kinetic network contributes net stoichiometric coefficients nu[A,i],
net progress rates omega[i], and participation flags delta[B,i]. The coupling
of A with B is |sum_i nu[A,i] omega[i] delta[B,i]| / max(P_A,C_A).
P_A and C_A are gross formation/consumption defined with signed net fluxes.
We propagate the maximum product of directed edge weights from each target.

This is standard CHEMICAL DRGEP, not the specialized P-DRGEP energy-target
procedure of Bellemans et al. (2021). See references in the README.
"""
import heapq
import numpy as np


def direct_interaction_coefficients(nu, net_rates, participation, *, eps=1e-300):
    """Return directed (A,B) interaction matrix with values in [0,1].

    nu: (species, reactions) product-minus-reactant stoichiometric matrix.
    net_rates: (reactions,) signed net progress rates in consistent units.
    participation: (species, reactions) bool; both reactants/products counted.
    Species with no production or consumption have a zero row.
    """
    nu = np.asarray(nu, dtype=float)
    w = np.asarray(net_rates, dtype=float)
    participation = np.asarray(participation, dtype=bool)
    if nu.ndim != 2 or w.shape != (nu.shape[1],) or participation.shape != nu.shape:
        raise ValueError('Expected nu and participation (species,reactions); rates (reactions,)')
    if not np.all(np.isfinite(nu)) or not np.all(np.isfinite(w)):
        raise ValueError('Non-finite stoichiometry or rates')
    contributions = nu * w[None, :]
    prod = np.maximum(contributions, 0).sum(axis=1)
    cons = np.maximum(-contributions, 0).sum(axis=1)
    denominator = np.maximum(prod, cons)
    # Numerator is signed sum BEFORE taking abs: cancellation is important.
    numerators = np.abs(contributions @ participation.astype(float).T)
    result = np.divide(numerators, denominator[:, None],
                       out=np.zeros_like(numerators), where=denominator[:, None]>eps)
    np.clip(result, 0., 1., out=result)
    np.fill_diagonal(result, 1.)
    return result


def propagate_importance(direct, targets):
    """Maximum-product path from any target to each species (0..1)."""
    direct = np.asarray(direct, dtype=float)
    n = direct.shape[0] if direct.ndim == 2 else 0
    if direct.shape != (n,n) or not np.all(np.isfinite(direct)) or np.any(direct<0) or np.any(direct>1+1e-12):
        raise ValueError('Expected square finite [0,1] direct-interaction matrix')
    targets = sorted(set(int(t) for t in targets))
    if not targets or any(t<0 or t>=n for t in targets):
        raise ValueError('At least one valid target index required')
    importance = np.zeros(n)
    heap=[]
    for t in targets:
        importance[t]=1.
        heapq.heappush(heap,(-1.,t))
    while heap:
        negative,a=heapq.heappop(heap)
        weight=-negative
        if weight<importance[a]-1e-14:
            continue
        for b in np.flatnonzero(direct[a] > 0):
            val=weight*direct[a,b]
            if val>importance[b]+1e-14:
                importance[b]=val
                heapq.heappush(heap,(-val,int(b)))
    return importance


def multi_state_importance(nu, rates_by_state, participation, targets):
    """Union of maximal importance over a reference ensemble of states."""
    arr=np.asarray(rates_by_state,dtype=float)
    if arr.ndim!=2 or arr.shape[0]==0:
        raise ValueError('rates_by_state must be (samples, reactions) nonempty')
    return np.maximum.reduce([propagate_importance(
        direct_interaction_coefficients(nu,w,participation),targets) for w in arr])


def choose_species(importance, names, threshold=0.05, required=()):
    """Return names passing the importance threshold plus mandatory species."""
    names=list(names)
    imp=np.asarray(importance)
    if imp.shape != (len(names),) or not 0<=threshold<=1:
        raise ValueError('Invalid importance or threshold')
    return [name for name,value in zip(names,imp) if value>=threshold or name in required]


def gross_flux_interaction_coefficients(nu, forward_rates, reverse_rates, participation, *, eps=1e-300):
    """A *diagnostic* graph based on gross forward/reverse flux magnitudes.

    Unlike standard net-rate DRGEP, this conservative screening avoids
    cancellation near equilibrium. It is NOT a replacement for calibrated
    mechanism reduction or for plasma P-DRGEP; dense networks may retain all.
    """
    nu=np.asarray(nu,float)
    f=np.asarray(forward_rates,float)
    r=np.asarray(reverse_rates,float)
    p=np.asarray(participation,bool)
    if (nu.ndim!=2 or f.shape!=(nu.shape[1],) or r.shape!=f.shape or
        p.shape!=nu.shape or not all(np.all(np.isfinite(a)) for a in (nu,f,r))
        or np.any(f<0) or np.any(r<0)):
        raise ValueError('Invalid gross-flux interaction arrays')
    magnitude=np.abs(nu)*(f+r)[None,:]
    denom=magnitude.sum(axis=1)
    numerator=magnitude@p.astype(float).T
    d=np.divide(numerator,denom[:,None],out=np.zeros_like(numerator),where=denom[:,None]>eps)
    np.clip(d,0.,1.,out=d)
    np.fill_diagonal(d,1.)
    return d

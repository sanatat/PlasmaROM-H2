"""Nonnegative data-fitted toy closure, with held-out evaluation.

Constrained least squares estimates a two-feature nonnegative effective
B=(F1+F2) -> R transfer rate from full surrogate trajectories. It is a deliberately
simple data-driven reduced model, not deep learning, PINN, or real PAC.
"""
import numpy as np
from scipy.optimize import nnls
from scipy.integrate import solve_ivp
from .kinetics import Parameters,simulate,pulse,qssa_partition


def fit_closure(train_conditions=((.5,.2),(.5,.55),(1.,.2),(1.,.55),(1.5,.2),(1.5,.55))):
    X=[];Y=[]
    for strength,width in train_conditions:
        p=Parameters(pulse_strength=strength,pulse_width=width)
        sol=simulate(p)
        drive=pulse(sol.t,p)
        # F2's exact output flux is deliberately fitted using full reference states.
        b=sol.y[1]+sol.y[2]
        X.append(np.column_stack((b*drive, b*drive**2)))
        Y.append(p.activated_oxidation*sol.y[2])
    coeff,residual=nnls(np.vstack(X),np.concatenate(Y))
    return coeff,float(residual)


def learned_rhs(t,y,p,coeff):
    f0,b,r,product=y
    drive=float(pulse(t,p))
    flux=max(0.,float(coeff[0]*drive+coeff[1]*drive*drive))*b
    q01=p.excitation*drive*f0
    # B is a physical F1+F2 sum; don't treat it as pure F1.
    q10=p.relaxation10*(1.-qssa_partition(t,p))*b
    q0r=p.background_oxidation*f0
    qrp=p.radical_oxidation*r
    qr0=p.radical_quench*r
    return np.array([-q01+q10-q0r+qr0, q01-q10-flux,
                     flux+q0r-qrp-qr0, qrp])


def simulate_learned(p,coeff,method='BDF'):
    coeff=np.asarray(coeff,float)
    if coeff.shape!=(2,) or not np.all(np.isfinite(coeff)) or np.any(coeff<0):
        raise ValueError('Expected 2 finite nonnegative closure weights')
    return solve_ivp(lambda t,y:learned_rhs(t,y,p,coeff),(0.,6.),[1.,0.,0.,0.],
                     t_eval=np.linspace(0.,6.,501),method=method,rtol=1e-8,atol=1e-11)


def held_out_metrics(coeff, cases=((2.5,.12),(6.,.34),(.8,.42))):
    rows=[]
    for strength,width in cases:
        p=Parameters(pulse_strength=strength,pulse_width=width)
        reference=simulate(p)
        qssa=simulate(p,reduced=True)
        pred=simulate_learned(p,coeff)
        if not (reference.success and pred.success and qssa.success):
            raise RuntimeError('Full, QSSA or learned integration failed')
        rows.append({'pulse_strength':strength,'pulse_width':width,
                     'learned_max_abs_product_error':float(np.max(np.abs(reference.y[4]-pred.y[3]))),
                     'qssa_max_abs_product_error':float(np.max(np.abs(reference.y[4]-qssa.y[3]))),
                     'learned_total_pool_residual':float(np.max(np.abs(pred.y.sum(axis=0)-1))),
                     'learned_min_pool':float(pred.y.min())})
    return rows

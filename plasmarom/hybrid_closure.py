"""Single-parameter bounded QSSA lag correction. Synthetic demo only.

Nonnegative lag time fitted on toy training trajectories. Relative to pure
QSSA, off-manifold behavior may worsen. Retains positive conservative rates.
"""
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.integrate import solve_ivp
from .kinetics import Parameters,pulse,rates,qssa_partition,simulate


def lag_fraction(t,p,tau):
    a01,a12=rates(t,p)
    # d/dt of second Gaussian pulse is analytic; avoid finite-difference bias.
    w=p.pulse_width
    deriv=p.pulse_strength*(
        -2*(t-1.2)/w**2*np.exp(-((t-1.2)/w)**2)
        -.75*2*(t-3.)/w**2*np.exp(-((t-3.)/w)**2))
    d_rate=p.second_excitation*deriv
    out=p.relaxation21+p.activated_oxidation
    denom=out+a12
    frac=float(qssa_partition(t,p))
    d_frac=float(out*d_rate/denom**2)
    return float(np.clip(frac-tau*d_frac,0,1))


def simulate_lag(p=Parameters(),tau=0.,n_points=201):
    if not np.isfinite(tau) or tau<0:raise ValueError('tau must be nonnegative')
    def rhs(t,y):
        f0,b,r,pr=y
        a01,_=rates(t,p)
        frac=lag_fraction(t,p,tau)
        f1=(1-frac)*b;f2=frac*b
        flux01=a01*f0;flux10=p.relaxation10*f1
        flux2r=p.activated_oxidation*f2;flux0r=p.background_oxidation*f0
        fluxrp=p.radical_oxidation*r;fluxr0=p.radical_quench*r
        return [-flux01+flux10-flux0r+fluxr0,flux01-flux10-flux2r,
                flux2r+flux0r-fluxrp-fluxr0,fluxrp]
    return solve_ivp(rhs,(0,6),[1.,0.,0.,0.],t_eval=np.linspace(0,6,n_points),
                     method='BDF',rtol=1e-8,atol=1e-11)


def fit_and_holdout(train=((.5,.2),(.5,.55),(1.5,.2),(1.5,.55)),
                    holdout=((2.5,.12),(6.,.34),(.8,.42)),
                    validation=((.8,.30),(2.,.45))):
    if (not train or not validation or not holdout or
        set(train)&set(validation) or set(train)&set(holdout) or set(validation)&set(holdout)):
        raise ValueError('Train, validation, holdout must be nonempty, disjoint conditions')
    references=[simulate(Parameters(pulse_strength=s,pulse_width=w),n_points=201)
                for s,w in train]
    def objective(tau):
        err=0
        for (s,w),ref in zip(train,references):
            fitted=simulate_lag(Parameters(pulse_strength=s,pulse_width=w),tau)
            err+=float(np.mean((ref.y[-1]-fitted.y[-1])**2))
        return err
    opt=minimize_scalar(objective,bounds=(0,.2),method='bounded',options={'xatol':1e-3})
    tau=float(opt.x)
    validation_rows=[]
    for s,w in validation:
        p=Parameters(pulse_strength=s,pulse_width=w)
        ref=simulate(p,n_points=201);qssa=simulate(p,reduced=True,n_points=201)
        lag=simulate_lag(p,tau)
        validation_rows.append({'pulse_strength':s,'pulse_width':w,
                                'lag_error':float(np.max(abs(ref.y[-1]-lag.y[-1]))),
                                'qssa_error':float(np.max(abs(ref.y[-1]-qssa.y[-1])))})
    # The fallback is selected ONLY from validation, never from the holdout.
    improved=max(z['lag_error'] for z in validation_rows) < max(z['qssa_error'] for z in validation_rows)
    selected='lag_fit' if improved else 'qssa_fallback'
    rows=[]
    for s,w in holdout:
        p=Parameters(pulse_strength=s,pulse_width=w)
        ref=simulate(p,n_points=201);qssa=simulate(p,reduced=True,n_points=201)
        lag=simulate_lag(p,tau)
        rows.append({'pulse_strength':s,'pulse_width':w,
                     'lag_max_abs_product_error':float(max(abs(ref.y[-1]-lag.y[-1]))),
                     'qssa_max_abs_product_error':float(max(abs(ref.y[-1]-qssa.y[-1]))),
                     'selected_max_abs_product_error':float(max(abs(ref.y[-1]-(lag.y[-1] if improved else qssa.y[-1])))),
                     'selected_model':selected,
                     'lag_mass_drift':float(max(abs(lag.y.sum(axis=0)-1))),
                     'lag_min_population':float(lag.y.min())})
    return {'fitted_tau_dimensionless':tau,'train_objective_mse_sum':float(opt.fun),
            'validation':validation_rows,'selected_model_by_validation':selected,
            'heldout':rows,'limitation':'Synthetic QSSA correction; no calibrated plasma and no claim of improved generalization'}

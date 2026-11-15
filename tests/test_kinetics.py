import numpy as np
from plasmarom.kinetics import Parameters, simulate, rhs_full, rhs_reduced

def test_flux_is_conservative():
    y=np.array([.7,.1,.02,.08,.1])
    assert abs(rhs_full(.9,y).sum())<1e-12
    assert abs(rhs_reduced(.9,np.array([.7,.1,.1,.1])).sum())<1e-12

def test_positive_and_conservation():
    for reduced in (False,True):
        a=simulate(reduced=reduced)
        assert a.success
        assert np.max(np.abs(a.y.sum(axis=0)-1))<1e-8
        assert a.y.min()>-1e-9

def test_reduced_output_at_baseline():
    a,b=simulate(),simulate(reduced=True)
    assert np.max(abs(a.y[4]-b.y[3]))<0.05
    assert np.max(abs(a.y[3]-b.y[2]))<0.05

def test_rk45_matches_bdf():
    a,b=simulate(method='BDF'),simulate(method='RK45')
    assert np.max(abs(a.y-b.y))<1e-6


def test_qssa_reconstruction_is_conservative_and_fast_pool_nonzero():
    from plasmarom.kinetics import qssa_reconstruct
    p=Parameters(pulse_strength=2.)
    red=simulate(p,reduced=True)
    rec=qssa_reconstruct(red,p)
    assert np.max(np.abs(rec.sum(axis=0)-1))<1e-8
    assert rec[2].max()>0


def test_radau_bdf_agree_in_strong_pulse():
    p=Parameters(pulse_strength=25.,pulse_width=.12)
    a=simulate(p,method='BDF')
    b=simulate(p,method='Radau')
    assert a.success and b.success
    assert np.max(np.abs(a.y-b.y))<1e-5


def test_lumped_qssa_rhs_agrees_with_full_fast_manifold():
    from plasmarom.kinetics import qssa_partition,full_observables,qssa_reconstruct
    p=Parameters(pulse_strength=5.,pulse_width=.12)
    for t in (.3,1.2,1.27,3.):
        B=.32;f0=.45;r=.12;prod=.11
        frac=qssa_partition(t,p)
        full=np.array([f0,B*(1-frac),B*frac,r,prod])
        dy=rhs_full(t,full,p)
        projected=np.array([dy[0],dy[1]+dy[2],dy[3],dy[4]])
        reduced_rhs=rhs_reduced(t,np.array([f0,B,r,prod]),p)
        np.testing.assert_allclose(reduced_rhs,projected,rtol=1e-12,atol=1e-12)
    sol=simulate(p,reduced=True)
    rec=qssa_reconstruct(sol,p)
    np.testing.assert_allclose(rec[1]+rec[2],sol.y[1],rtol=1e-12,atol=1e-12)
    full=simulate(p)
    obs=full_observables(full)
    assert np.max(abs(obs[1]-sol.y[1]))<.04


def test_strong_pulse_product_error_regression():
    p=Parameters(pulse_strength=100.,pulse_width=.12)
    a=simulate(p);b=simulate(p,reduced=True)
    assert a.success and b.success
    assert np.max(np.abs(a.y[4]-b.y[3]))<.08

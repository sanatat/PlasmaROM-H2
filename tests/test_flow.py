import numpy as np
from plasmarom.flow import SpectralFlow, FlowConfig, run_flow

def test_projection():
    f=SpectralFlow(24)
    f.u+=0.02*np.cos(2*f.x)
    f.u,f.v=f.project(f.u,f.v)
    assert np.max(abs(f.divergence()))<1e-10

def test_taylor_green_exact():
    f=SpectralFlow(24,nu=0.01)
    for _ in range(20):f.step(0.02)
    assert f.taylor_green_error(0.4)<1e-3

def test_reactive_advection_sanity():
    d=run_flow(FlowConfig(n=20,dt=0.02,t_end=0.6))
    assert d['peak_divergence']<1e-9
    assert d['sum_deviation']<1e-8
    assert d['min_pool']>-1e-8


def test_strong_pulse_stays_nonnegative_and_conservative():
    from plasmarom.kinetics import Parameters
    p=Parameters(pulse_strength=100.)
    for reduced in (False,True):
        d=run_flow(FlowConfig(n=20,dt=0.02,t_end=1.5,reduced=reduced),p)
        assert d['min_pool']>=-1e-10
        assert d['sum_deviation']<1e-10
        assert d['max_chemical_substeps_per_half_step']>1


def test_reject_transport_cfl_and_bad_params():
    import pytest
    with pytest.raises(ValueError,match='CFL'):
        run_flow(FlowConfig(n=20,dt=0.5,t_end=0.5))
    with pytest.raises(ValueError):
        run_flow(FlowConfig(n=15))


def test_refinement_agreement_2d():
    coarse=run_flow(FlowConfig(n=20,dt=0.03,t_end=0.6))
    fine=run_flow(FlowConfig(n=20,dt=0.015,t_end=0.6))
    assert np.mean(abs(coarse['species'][-1]-fine['species'][-1]))<0.005


def test_negative_chemistry_coefficients_are_rejected():
    import pytest
    from plasmarom.kinetics import Parameters
    with pytest.raises(ValueError,match='chemical transfer'):
        run_flow(FlowConfig(n=20,dt=0.02,t_end=.1),Parameters(activated_oxidation=-1.))


def test_finite_volume_individual_scalar_conservation():
    """Regression: v0.2 lost up to 0.24% of individual species at 64²."""
    from plasmarom.flow import (divergence_free_face_velocities,
        scalar_transport,finite_volume_cfl)
    for n in (16,32,64):
        flow=SpectralFlow(n)
        x,y=flow.x,flow.y
        gaussian=np.exp(-((x-2.)**2+(y-3.)**2)/.4)
        scalars=np.stack([gaussian,1-gaussian])
        before=scalars.sum(axis=(1,2)).copy()
        faces=divergence_free_face_velocities(flow.u,flow.v,flow.dx)
        assert .005*finite_volume_cfl(faces,flow.dx,0.) < .9
        fvdiv=(faces[0]-np.roll(faces[0],1,axis=0)+
               faces[1]-np.roll(faces[1],1,axis=1))/flow.dx
        assert np.max(abs(fvdiv))<1e-11
        for step in range(100):
            scalars += .005*scalar_transport(scalars,flow.u,flow.v,flow.dx,0.,faces)
        assert np.max(abs(scalars.sum(axis=(1,2))-before))<1e-9
        assert scalars.min()>-1e-14
        assert np.max(abs(scalars.sum(axis=0)-1))<1e-12


def test_constant_scalar_advects_with_zero_change():
    from plasmarom.flow import scalar_transport
    flow=SpectralFlow(32)
    phi=np.ones((3,32,32))*np.array([.2,.3,.5])[:,None,None]
    assert np.max(abs(scalar_transport(phi,flow.u,flow.v,flow.dx,.01)))<1e-12

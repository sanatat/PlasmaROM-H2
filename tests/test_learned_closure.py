import numpy as np
from plasmarom.kinetics import Parameters
from plasmarom.learned_closure import fit_closure,simulate_learned,held_out_metrics


def test_nonnegative_constrained_learned_closure():
    w,residual=fit_closure()
    assert w.shape==(2,) and np.all(w>=0)
    assert residual>=0
    sol=simulate_learned(Parameters(pulse_strength=3.),w)
    assert sol.success
    assert sol.y.min()>=-1e-9
    assert np.max(np.abs(sol.y.sum(axis=0)-1))<1e-8


def test_out_of_sample_closure_results_are_finite():
    w,_=fit_closure()
    for row in held_out_metrics(w):
        assert np.isfinite(row['learned_max_abs_product_error'])
        assert row['learned_total_pool_residual']<1e-7

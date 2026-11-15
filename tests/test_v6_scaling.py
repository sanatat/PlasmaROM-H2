import pytest
from plasmarom.parallel import fixed_workload_scaling

pytestmark = pytest.mark.slow


def test_repeated_scaling_is_numerically_reproducible(tmp_path):
    d=fixed_workload_scaling(tmp_path,powers=(100.,500.),workers=(1,2),repetitions=2)
    assert len(d['by_workers']['1']['wall_times_s'])==2
    assert len(d['by_workers']['2']['wall_times_s'])==2
    assert d['by_workers']['1']['speedup_vs_serial']==pytest.approx(1.)
    assert d['by_workers']['2']['median_wall_s']>0
    assert (tmp_path/'fixed_workload_scaling.json').exists()


def test_invalid_scaling_parameters(tmp_path):
    for kw in ({'workers':(2,)},{'repetitions':0},{'workers':(1,1)}):
        with pytest.raises(ValueError):
            fixed_workload_scaling(tmp_path,**kw)

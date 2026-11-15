import json
import pytest
from plasmarom import cantera_flame


def test_grid_study_passes_only_small_speed_change(monkeypatch,tmp_path):
    def simulate(out,**kwargs):
        slope=kwargs['refine_slope']
        return {'flame_speed_m_s':2.00 if slope>.1 else 2.03,
                'grid_points':50 if slope>.1 else 150}
    monkeypatch.setattr(cantera_flame,'run_flame',simulate)
    result=cantera_flame.run_flame_refinement_study(tmp_path,tolerance=.02)
    assert result['status']=='GRID_SENSITIVITY_PASS'
    assert result['relative_flame_speed_change']==pytest.approx(.03/2.03)
    assert (tmp_path/'grid_refinement_summary.json').exists()


def test_grid_study_reports_failure_and_no_fake_success(monkeypatch,tmp_path):
    def simulate(out,**kwargs):
        return {'flame_speed_m_s':2 if kwargs['refine_slope']>.1 else 3,
                'grid_points':30}
    monkeypatch.setattr(cantera_flame,'run_flame',simulate)
    result=cantera_flame.run_flame_refinement_study(tmp_path,tolerance=.05)
    assert result['status']=='GRID_SENSITIVITY_FAIL'
    with pytest.raises(ValueError):
        cantera_flame.run_flame_refinement_study(tmp_path,tolerance=0)

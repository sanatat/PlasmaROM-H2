import pytest
from plasmarom.hybrid_closure import fit_and_holdout


def test_validation_only_selects_model_without_holdout_leakage():
    data=fit_and_holdout()
    assert data['selected_model_by_validation'] in {'lag_fit','qssa_fallback'}
    for row in data['heldout']:
        if data['selected_model_by_validation']=='lag_fit':
            assert row['selected_max_abs_product_error']==row['lag_max_abs_product_error']
        else:
            assert row['selected_max_abs_product_error']==row['qssa_max_abs_product_error']


def test_data_split_overlap_is_rejected():
    with pytest.raises(ValueError,match='disjoint'):
        fit_and_holdout(validation=((.5,.2),))

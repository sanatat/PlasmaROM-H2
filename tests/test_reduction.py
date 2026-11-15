import numpy as np
import pytest
from plasmarom.reduction import (direct_interaction_coefficients,
    propagate_importance,multi_state_importance,choose_species)


def test_direct_interaction_definition_and_cancellation():
    # A produced in reaction 1 and consumed by reaction 2, B participates in
    # both: numerator cancels despite nonzero gross production/consumption.
    nu=np.array([[1,-1],[0,0],[-1,1]])
    participation=np.array([[1,1],[1,1],[1,1]],dtype=bool)
    c=direct_interaction_coefficients(nu,np.array([1.,1.]),participation)
    assert c[0,1]==pytest.approx(0.)
    assert c[0,0]==1.


def test_max_product_path_and_target_union():
    g=np.eye(4)
    g[0,1]=0.8;g[1,2]=0.5;g[0,2]=0.1;g[2,3]=0.7
    importance=propagate_importance(g,[0])
    assert importance==pytest.approx([1.,0.8,0.4,0.28])
    new=propagate_importance(g,[0,3])
    assert new[3]==pytest.approx(1.)


def test_select_and_multistate():
    nu=np.array([[-1,0],[1,-1],[0,1]],float)
    participation=np.abs(nu)>0
    rates=np.array([[1.,0.],[0.,1.]])
    imp=multi_state_importance(nu,rates,participation,targets=[0])
    assert len(imp)==3
    assert choose_species([1.,.4,.02],['A','B','C'],.2,required=('C',))==['A','B','C']


def test_shapes_and_invalid_inputs():
    with pytest.raises(ValueError):
        direct_interaction_coefficients(np.ones((3,2)),np.ones(3),np.ones((3,2)))
    with pytest.raises(ValueError):
        propagate_importance(np.ones((2,2)),targets=[])

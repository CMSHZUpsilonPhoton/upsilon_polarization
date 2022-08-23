import pytest
import numpy as np

from upsilon_polarization import full_method

cms_data = full_method.load_CMS_data()

def test_get_weight_full_method_scalars():
    assert full_method.get_weight(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1) == -4.071550835635227

a = np.ones(60)
def test_get_weight_full_method_arrays():
    assert np.array_equal(full_method.get_weight(a, a, a, a, a, a, a, a, a, a, a), np.ones(60)*(full_method.get_weight(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)))

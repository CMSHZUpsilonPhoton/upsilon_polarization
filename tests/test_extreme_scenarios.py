import pytest
import math

from upsilon_polarization import extreme_scenarios

import ROOT

def test_upsilon_polarization_cos_angle():
    assert extreme_scenarios.get_cos_angle(1, 2, 3, 1, 1, 3) == 0

def test_get_weight_extreme_scenarios_good_cos_angle():
     assert extreme_scenarios.get_weight(0,1) == 0.75


def test_get_weight_extreme_scenarios_bad_cos_angle():
    with pytest.raises(ValueError) as exc_info:
        extreme_scenarios.get_weight(3, -1)

def test_get_weight_extreme_scenarios_bad_syst():
    with pytest.raises(ValueError) as exc_info:
        extreme_scenarios.get_weight(3, -10)

import pytest

from upsilon_polarization import full_method

cms_data = full_method.load_CMS_data()

def test_get_weight_extreme_scenarios_good_cos_angle():
    assert full_method.get_weight(1, 1, 1, 1, 1, 1, 1, 1, 1, cms_data, 1, 1, 0) == -4.073500193034304


def test_get_weight_full_method_bad_syst():
    with pytest.raises(ValueError) as exc_info:
        full_method.get_weight(1, 1, 1, 1, 1, 1, 1, 1, 1, cms_data, 1, 1, -10)

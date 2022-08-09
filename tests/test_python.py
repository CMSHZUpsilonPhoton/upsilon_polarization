from ctypes import CDLL, c_double, c_int

upslib = CDLL("lib/upsilon_polarization.so")

upslib.upsilon_polarization_cos_angle.argtypes = [c_double]
upslib.upsilon_polarization_cos_angle.restype = c_double

print(upslib.upsilon_polarization_cos_angle(123.0))
print(f"type: {type(upslib.upsilon_polarization_cos_angle(123.0))}")


upslib.get_weight_extreme_scenarios.argtypes = [c_double, c_int]
upslib.get_weight_extreme_scenarios.restype = c_double

print(upslib.get_weight_extreme_scenarios(123.0, -1))
print(f"type: {type(upslib.get_weight_extreme_scenarios(123.0, -1))}")



def run_upsilon_polarization_cos_angle(pt):
    return upslib.upsilon_polarization_cos_angle(pt)

def run_get_weight_extreme_scenarios(cos_angle, syst):
    return upslib.get_weight_extreme_scenarios(cos_angle, syst)




def test_upsilon_polarization_cos_angle():
    assert run_upsilon_polarization_cos_angle(3) == 3

    
def test_get_weight_extreme_scenarios():
    assert run_get_weight_extreme_scenarios(3, -1) == -12
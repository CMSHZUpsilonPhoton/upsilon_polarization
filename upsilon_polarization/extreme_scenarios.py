from ctypes import CDLL, c_double, c_int, c_char_p

upslib = CDLL("upsilon_polarization/lib/upsilon_polarization.so")

upslib.upsilon_polarization_cos_angle.argtypes = [c_double, c_double, c_double, c_double, c_double, c_double]
upslib.upsilon_polarization_cos_angle.restype = c_double


upslib.get_weight_extreme_scenarios.argtypes = [c_double, c_int]
upslib.get_weight_extreme_scenarios.restype = c_double

__all__ = ["get_weight"]

def get_cos_angle(pt_plus, eta_plus, phi_plus, pt_minus, eta_minus, phi_minus):
    return upslib.upsilon_polarization_cos_angle(pt_plus, eta_plus, phi_plus, pt_minus, eta_minus, phi_minus)

def get_weight(cos_angle: float, syst: int) -> float:
    if syst != -1 and syst != 0 and syst != +1:
        raise ValueError(f"The provided value for 'syst' ({syst}) is not valid. It should be -1, 0 or +1.")

    if cos_angle > +1 or cos_angle < -1:
        raise ValueError(f"The provided value for 'cos_angle' ({cos_angle}) is not valid. It should between -1, and +1.")
        
    return upslib.get_weight_extreme_scenarios(cos_angle, syst)
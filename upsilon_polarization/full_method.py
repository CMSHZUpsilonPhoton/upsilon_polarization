from ctypes import CDLL, c_double, c_int, c_char_p, Structure, c_void_p, byref, POINTER

upslib = CDLL("upsilon_polarization/lib/upsilon_polarization.so")


class CoefficientsHists(Structure):
    _fields_ = [
        ("A0_low_y", c_void_p),
        ("A1_low_y", c_void_p),
        ("A2_low_y", c_void_p),
        ("A0_high_y", c_void_p),
        ("A1_high_y", c_void_p),
        ("A2_high_y", c_void_p),
        ("A0_low_y_stat", c_void_p),
        ("A1_low_y_stat", c_void_p),
        ("A2_low_y_stat", c_void_p),
        ("A0_high_y_stat", c_void_p),
        ("A1_high_y_stat", c_void_p),
        ("A2_high_y_stat", c_void_p),
        ("A0_low_y_syst", c_void_p),
        ("A1_low_y_syst", c_void_p),
        ("A2_low_y_syst", c_void_p),
        ("A0_high_y_syst", c_void_p),
        ("A1_high_y_syst", c_void_p),
        ("A2_high_y_syst", c_void_p),
    ]


upslib.load_CMS_data.argtypes = []
upslib.load_CMS_data.restype = CoefficientsHists




upslib.get_weight_full_method.argtypes = [
    c_double,
    c_double,
    c_double,
    c_double,
    c_double,
    c_double,
    c_double,
    c_double,
    c_double,
    POINTER(CoefficientsHists),
    c_double,
    c_double,
    c_int,
]
upslib.get_weight_full_method.restype = c_double



__all__ = ["load_CMS_data", "get_weight_full_method"]


def load_CMS_data():
    cms_data = upslib.load_CMS_data()
    return cms_data


def get_weight(
    pt_mu_plus,
    eta_mu_plus,
    phi_mu_plus,
    pt_mu_minus,
    eta_mu_minus,
    phi_mu_minus,
    pt_photon,
    eta_photon,
    phi_photon,
    coefficients_hists,
    z_rapidity,
    qT,
    syst,
) -> float:
    if syst != -1 and syst != 0 and syst != +1:
        raise ValueError(
            f"The provided value for 'syst' ({syst}) is not valid. It should be -1 (DOWN), 0 (NOMINAL) or +1 (UP)."
        )

    weight = upslib.get_weight_full_method(
        pt_mu_plus,
        eta_mu_plus,
        phi_mu_plus,
        pt_mu_minus,
        eta_mu_minus,
        phi_mu_minus,
        pt_photon,
        eta_photon,
        phi_photon,
        byref(coefficients_hists),
        z_rapidity,
        qT,
        syst,
    )

    
    return weight
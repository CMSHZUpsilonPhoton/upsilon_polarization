from ctypes import CDLL, c_double, c_int, c_void_p

from numba import vectorize, double, int32, float32, float64


upslib = CDLL("upsilon_polarization/lib/upsilon_polarization.so")

upslib.load_CMS_data.argtypes = []
upslib.load_CMS_data.restype = c_void_p


upslib_get_weight_full_method = upslib.get_weight_full_method
upslib_get_weight_full_method.argtypes = [
    c_double,
    c_double,
    c_double,
    c_double,
    c_double,
    c_double,
    c_double,
    c_double,
    c_double,
    c_void_p,
    c_double,
    c_double,
    c_int,
]
upslib_get_weight_full_method.restype = c_double


def load_CMS_data():
    cms_data = upslib.load_CMS_data()
    return cms_data

_cms_data = load_CMS_data()


@vectorize([float64(float64, float64, float64, float64, float64, float64, float64, float64, float64, float64, float64)])
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
    z_rapidity,
    qT) -> float:

    weight = upslib_get_weight_full_method(
        pt_mu_plus,
        eta_mu_plus,
        phi_mu_plus,
        pt_mu_minus,
        eta_mu_minus,
        phi_mu_minus,
        pt_photon,
        eta_photon,
        phi_photon,
        _cms_data,
        z_rapidity,
        qT,
        0,
    )

    return weight

@vectorize([float64(float64, float64, float64, float64, float64, float64, float64, float64, float64, float64, float64)])
def get_weight_plus(
    pt_mu_plus,
    eta_mu_plus,
    phi_mu_plus,
    pt_mu_minus,
    eta_mu_minus,
    phi_mu_minus,
    pt_photon,
    eta_photon,
    phi_photon,
    z_rapidity,
    qT) -> float:


    weight = upslib_get_weight_full_method(
        pt_mu_plus,
        eta_mu_plus,
        phi_mu_plus,
        pt_mu_minus,
        eta_mu_minus,
        phi_mu_minus,
        pt_photon,
        eta_photon,
        phi_photon,
        _cms_data,
        z_rapidity,
        qT,
        1,
    )

    return weight

@vectorize([float64(float64, float64, float64, float64, float64, float64, float64, float64, float64, float64, float64)])
def get_weight_minus(
    pt_mu_plus,
    eta_mu_plus,
    phi_mu_plus,
    pt_mu_minus,
    eta_mu_minus,
    phi_mu_minus,
    pt_photon,
    eta_photon,
    phi_photon,
    z_rapidity,
    qT) -> float:


    weight = upslib_get_weight_full_method(
        pt_mu_plus,
        eta_mu_plus,
        phi_mu_plus,
        pt_mu_minus,
        eta_mu_minus,
        phi_mu_minus,
        pt_photon,
        eta_photon,
        phi_photon,
        _cms_data,
        z_rapidity,
        qT,
        -1,
    )

    return weight


__all__ = ["load_CMS_data", "get_weight", "get_weight_plus", "get_weight_minus"]

#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "Math/VectorUtil.h"

extern "C"
{
    double upsilon_polarization_angle(double pt_plus, double eta_plus, double phi_plus, double pt_minus, double eta_minus, double phi_minus);
    double get_weight_extreme_scenarios(double cos_angle, int syst);
}
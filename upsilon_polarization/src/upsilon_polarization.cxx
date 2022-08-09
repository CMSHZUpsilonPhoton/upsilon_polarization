#include "upsilon_polarization.h"

using namespace ROOT::Math;

extern "C"
{
	const double MUON_MASS = 0.1056583755;

	double upsilon_polarization_cos_angle(double pt_plus, double eta_plus, double phi_plus, double pt_minus, double eta_minus, double phi_minus)
	{
		auto gen_mu_plus = PtEtaPhiMVector(pt_plus, eta_plus, phi_plus, MUON_MASS);
		auto gen_mu_minus = PtEtaPhiMVector(pt_minus, eta_minus, phi_minus, MUON_MASS);

		auto upsilon_CM = (gen_mu_plus + gen_mu_minus);

		auto gen_mu_plus_upsilon_frame = VectorUtil::boost(gen_mu_plus, -upsilon_CM);

		return VectorUtil::CosTheta(gen_mu_plus_upsilon_frame, upsilon_CM);
	}

	double get_weight_extreme_scenarios(double cos_angle, int syst)
	{
		// check syst input
		if (syst != -1 && syst != 0 && syst != +1)
		{
			printf("[ERROR] [POLARIZATION WEIGHT] Invalid argument. \"syst\" should be -1, 0 (unpolarized) or +1.");
			// throw std::invalid_argument("ERROR");
			exit(-1);
		}

		// polarization weight
		auto polWgt_Nominal = 1.0;								 // unpolarized
		auto polWgt_PLUS = (3. / 4.) * (1. + pow(cos_angle, 2));	 // transverse
		auto polWgt_MINUS = (3. / 2.) * (1. - pow(cos_angle, 2)); // longitudinal}

		double polarization_weight = 0.;
		if (syst == -1)
		{
			polarization_weight = polWgt_MINUS;
		}
		else if (syst == +1)
		{
			polarization_weight = polWgt_PLUS;
		}
		else
		{
			polarization_weight = polWgt_Nominal;
		}

		return polarization_weight;
	}
}
#include "upsilon_polarization.h"

using namespace ROOT::Math;

extern "C"
{
	double upsilon_polarization_cos_angle(double pt)
	{
		// TVector3 boostVector = -1*(genMuPlus + genMuMinus).BoostVector();
		// TLorentzVector genMuPlusUpsilonFrame(genMuPlus);
		// genMuPlusUpsilonFrame.Boost(boostVector);
		// return genMuPlusUpsilonFrame.Vect().Unit().Dot((genMuPlus + genMuMinus).Vect().Unit());
		auto b = PtEtaPhiMVector(pt, 2, 3, 4);
		return b.pt();
	}

	double get_weight_extreme_scenarios(double cos_angle, int syst)
	{
		// check syst input
		if (syst != -1 && syst != 0 && syst != +1)
		{
			throw std::invalid_argument("[ERROR] Invalid argument. \"syst\" should be -1, 0 (unpolarized) or +1.");
		}

		// polarization weight
		auto cosAngle = upsilon_polarization_cos_angle(cos_angle);
		auto polWgt_Nominal = 1.0;								// unpolarized
		auto polWgt_PLUS = (3. / 4.) * (1. + pow(cosAngle, 2));	// transverse
		auto polWgt_MINUS = (3. / 2.) * (1. - pow(cosAngle, 2)); // longitudinal}

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
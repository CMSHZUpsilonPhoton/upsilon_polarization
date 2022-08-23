// #include "upsilon_polarization.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "Math/VectorUtil.h"

#include "TMath.h"
#include "TH1F.h"
#include "TMemFile.h"

#include "cms_Z_polarization_data.hpp"

using namespace ROOT::Math;

extern "C"
{
	constexpr double MUON_MASS = 0.1056583755;
	constexpr double PROTON_MASS = 0.93827208816;
	constexpr double PROTON_ENERGY = 13000.0;

	/////////////////////////////////////////////////////
	/////////////////////////////////////////////////////
	// EXTREME
	/////////////////////////////////////////////////////
	/////////////////////////////////////////////////////

	double upsilon_polarization_cos_angle(double pt_mu_plus, double eta_mu_plus, double phi_mu_plus, double pt_mu_minus, double eta_mu_minus, double phi_mu_minus)
	{
		const auto gen_mu_plus = PtEtaPhiMVector(pt_mu_plus, eta_mu_plus, phi_mu_plus, MUON_MASS);
		const auto gen_mu_minus = PtEtaPhiMVector(pt_mu_minus, eta_mu_minus, phi_mu_minus, MUON_MASS);

		const auto upsilon_CM = (gen_mu_plus + gen_mu_minus);

		const auto gen_mu_plus_upsilon_frame = VectorUtil::boost(gen_mu_plus, upsilon_CM.BoostToCM());

		return VectorUtil::CosTheta(gen_mu_plus_upsilon_frame, upsilon_CM);
	}

	double get_weight_extreme_scenarios(double cos_angle, int syst)
	{
		// check syst input
		if (syst != -1 && syst != 0 && syst != +1)
		{
			printf("[ERROR] [POLARIZATION WEIGHT] Invalid argument. \"syst\" should be -1 (longitudinal), 0 (unpolarized) or +1 (transverse).");
			exit(-1);
		}

		// polarization weight
		const auto polWgt_Nominal = 1.0;								// unpolarized
		const auto polWgt_PLUS = (3. / 4.) * (1. + pow(cos_angle, 2));	// transverse
		const auto polWgt_MINUS = (3. / 2.) * (1. - pow(cos_angle, 2)); // longitudinal

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

	/////////////////////////////////////////////////////
	/////////////////////////////////////////////////////
	// ALTERNATIVE
	/////////////////////////////////////////////////////
	/////////////////////////////////////////////////////

	double get_cos_angle_alternative_method(
		const double &pt_mu_plus,
		const double &eta_mu_plus,
		const double &phi_mu_plus,
		const double &pt_mu_minus,
		const double &eta_mu_minus,
		const double &phi_mu_minus,
		const double &pt_photon,
		const double &eta_photon,
		const double &phi_photon)
	{
		// define particles
		const auto gen_mu_plus = PtEtaPhiMVector(pt_mu_plus, eta_mu_plus, phi_mu_plus, MUON_MASS);
		const auto gen_mu_minus = PtEtaPhiMVector(pt_mu_minus, eta_mu_minus, phi_mu_minus, MUON_MASS);
		const auto upsilon = (gen_mu_plus + gen_mu_minus);
		const auto photon = PtEtaPhiMVector(pt_photon, eta_photon, phi_photon, 0.0);
		const auto boson = (upsilon + photon);

		// Upsilon in the Z rest frame
		const auto upsilon_Z_rest_frame = VectorUtil::boost(upsilon, boson.BoostToCM());

		// positive muon in the Upsilon rest frame
		const auto gen_mu_plus_upsilon_rest_frame = VectorUtil::boost(gen_mu_plus, upsilon.BoostToCM());

		return VectorUtil::CosTheta(upsilon_Z_rest_frame, gen_mu_plus_upsilon_rest_frame);
	}

	double get_weight_alternative_method(
		const double &pt_mu_plus,
		const double &eta_mu_plus,
		const double &phi_mu_plus,
		const double &pt_mu_minus,
		const double &eta_mu_minus,
		const double &phi_mu_minus,
		const double &pt_photon,
		const double &eta_photon,
		const double &phi_photon,
		const int syst)
	{
		// check syst input
		if (syst != -1 && syst != 0 && syst != +1)
		{
			printf("[ERROR] [POLARIZATION WEIGHT] Invalid argument. \"syst\" should be -1 (longitudinal), 0 (unpolarized) or +1 (transverse).");
			exit(-1);
		}

		const double cos_angle = get_cos_angle_alternative_method(pt_mu_plus, eta_mu_plus, phi_mu_plus, pt_mu_minus, eta_mu_minus, phi_mu_minus, pt_photon, eta_photon, phi_photon);

		// polarization weight
		return 1.0 - (1.0 / 3.0) * cos_angle;
	}

	/////////////////////////////////////////////////////
	/////////////////////////////////////////////////////
	// FULL
	/////////////////////////////////////////////////////
	/////////////////////////////////////////////////////

	// Full Method
	// REF: https://indico.cern.ch/event/723934/contributions/2977415/attachments/1646018/2630726/PFaccioli_CMS_8May2018.pdf#search=Pietro%20Faccioli

	// math helper
	double pow2(double x) { return TMath::Power(x, 2); }
	double sqrt(double x) { return TMath::Sqrt(x); }
	double cos(double x) { return TMath::Cos(x); }
	double sin(double x) { return TMath::Sin(x); }
	double cos2(double x) { return pow2(cos(x)); }
	double sin2(double x) { return pow2(sin(x)); }
	double asin(double x) { return TMath::ASin(x); }
	double acos(double x) { return TMath::ACos(x); }

	// ROOT TH1 helper
	double get_bin_content(const double &value, TH1F *hist)
	{
		return hist->GetBinContent(hist->GetXaxis()->FindBin(value));
	}

	typedef struct CSAngles
	{
		double BIG_PHI;
		double cos_BIG_THETA;
		double BIG_THETA;
		double phi;
		double cos_theta;
		double theta;
	} CSAngles;

	CSAngles get_CSAngles(const double &pt_mu_plus,
						  const double &eta_mu_plus,
						  const double &phi_mu_plus,
						  const double &pt_mu_minus,
						  const double &eta_mu_minus,
						  const double &phi_mu_minus,
						  const double &pt_photon,
						  const double &eta_photon,
						  const double &phi_photon)
	{
		// define particles
		const auto gen_mu_plus = PtEtaPhiMVector(pt_mu_plus, eta_mu_plus, phi_mu_plus, MUON_MASS);
		const auto gen_mu_minus = PtEtaPhiMVector(pt_mu_minus, eta_mu_minus, phi_mu_minus, MUON_MASS);
		const auto upsilon = (gen_mu_plus + gen_mu_minus);
		const auto photon = PtEtaPhiMVector(pt_photon, eta_photon, phi_photon, 0.0);
		const auto boson = (upsilon + photon);

		// define beam beam and target
		const double PROTON_MOMEMTUM = sqrt(pow2(PROTON_ENERGY) - pow2(PROTON_MASS));
		const auto proton_plus = PxPyPzEVector(0., 0., PROTON_MOMEMTUM, PROTON_ENERGY);
		const auto proton_minus = PxPyPzEVector(0., 0., -PROTON_MOMEMTUM, PROTON_ENERGY);

		PxPyPzEVector beam;	  // beam
		PxPyPzEVector target; // target

		// looks like the CMS paper definition is swapped...
		if (boson.pz() >= 0)
		{
			beam = proton_plus;
			target = proton_minus;
		}
		else
		{
			beam = proton_minus;
			target = proton_plus;
		}

		// // just trying to swap the beam x target definition
		// if (boson.pz() < 0)
		// {
		// 	beam = proton_plus;
		// 	target = proton_minus;
		// }
		// else
		// {
		// 	beam = proton_minus;
		// 	target = proton_plus;
		// }

		// beta to boost to CS frame
		const auto beta_boson = boson.BoostToCM();

		// boost to CS frame
		const auto upsilon_CS_frame = VectorUtil::boost(upsilon, beta_boson);
		// auto photon_CS_frame = VectorUtil::boost(photon, beta_boson); // not used
		const auto beam_CS_frame = VectorUtil::boost(beam, beta_boson);
		const auto target_CS_frame = VectorUtil::boost(target, beta_boson);

		// Determine x, y and z axis for the CS frame
		const auto z_axis_CS_frame = (beam_CS_frame.Vect().Unit() - target_CS_frame.Vect().Unit()).Unit();
		const auto y_axis_CS_frame = (beam_CS_frame.Vect().Unit().Cross(-target_CS_frame.Vect().Unit())).Unit();
		const auto x_axis_CS_frame = y_axis_CS_frame.Cross(z_axis_CS_frame).Unit();

		// BIG_PHI
		// https://en.wikipedia.org/wiki/Atan2#/media/File:Arctangent2.svg
		const double BIG_PHI = TMath::ATan2((upsilon_CS_frame.Vect()).Dot(y_axis_CS_frame), ((upsilon_CS_frame.Vect()).Dot(x_axis_CS_frame)));

		// cos_BIG_THETA
		const double cos_BIG_THETA = z_axis_CS_frame.Dot((upsilon_CS_frame.Vect()).Unit());
		// Theta CS is not properly defined for Like-Sign (?)
		// if (cos_BIG_THETA < 0)
		// 	cos_BIG_THETA = -cos_BIG_THETA;
		const double BIG_THETA = acos(cos_BIG_THETA);

		// rigid translation to the upsilon CM frame
		// z axis is preserved
		const auto gen_mu_plus_upsilon_frame = VectorUtil::boost(gen_mu_plus, upsilon.BoostToCM());
		// auto gen_mu_minus_upsilon_frame = VectorUtil::boost(gen_mu_minus, upsilon.BoostToCM()); // not used

		// phi
		const double phi = TMath::ATan2((gen_mu_plus_upsilon_frame.Vect()).Dot(y_axis_CS_frame), ((gen_mu_plus_upsilon_frame.Vect()).Dot(x_axis_CS_frame)));

		// cos_theta
		const double cos_theta = z_axis_CS_frame.Dot((gen_mu_plus_upsilon_frame.Vect()).Unit());
		// Theta CS is not properly defined for Like-Sign  (?)
		// if (cos_theta < 0)
		// 	cos_theta = -cos_theta;
		const double theta = acos(cos_theta);

		auto cs_angles = CSAngles();
		cs_angles.BIG_PHI = BIG_PHI;
		cs_angles.cos_BIG_THETA = cos_BIG_THETA;
		cs_angles.BIG_THETA = BIG_THETA;
		cs_angles.phi = phi;
		cs_angles.cos_theta = cos_theta;
		cs_angles.theta = theta;

		return cs_angles;
	}

	typedef struct CoefficientsHists
	{
		void *A0_low_y;
		void *A1_low_y;
		void *A2_low_y;
		void *A0_high_y;
		void *A1_high_y;
		void *A2_high_y;
		void *A0_low_y_stat;
		void *A1_low_y_stat;
		void *A2_low_y_stat;
		void *A0_high_y_stat;
		void *A1_high_y_stat;
		void *A2_high_y_stat;
		void *A0_low_y_syst;
		void *A1_low_y_syst;
		void *A2_low_y_syst;
		void *A0_high_y_syst;
		void *A1_high_y_syst;
		void *A2_high_y_syst;
	} CoefficientsHists;

	void *load_CMS_data()
	{

		auto cms_Z_polarization_data_ = static_cast<char *>(static_cast<void *>(cms_Z_polarization_data));
		auto file = TMemFile("hepdata", cms_Z_polarization_data_, cms_Z_polarization_data_len);

		// nominal
		static auto A0_low_y_hist_ = (*file.Get<TH1F>("Table 1/Hist1D_y1"));
		static auto A1_low_y_hist_ = (*file.Get<TH1F>("Table 1/Hist1D_y2"));
		static auto A2_low_y_hist_ = (*file.Get<TH1F>("Table 1/Hist1D_y3"));
		static auto A0_high_y_hist_ = (*file.Get<TH1F>("Table 2/Hist1D_y1"));
		static auto A1_high_y_hist_ = (*file.Get<TH1F>("Table 2/Hist1D_y2"));
		static auto A2_high_y_hist_ = (*file.Get<TH1F>("Table 2/Hist1D_y3"));

		// stat
		static auto A0_low_y_stat_hist_ = (*file.Get<TH1F>("Table 1/Hist1D_y1_e1"));
		static auto A1_low_y_stat_hist_ = (*file.Get<TH1F>("Table 1/Hist1D_y2_e1"));
		static auto A2_low_y_stat_hist_ = (*file.Get<TH1F>("Table 1/Hist1D_y3_e1"));
		static auto A0_high_y_stat_hist_ = (*file.Get<TH1F>("Table 2/Hist1D_y1_e1"));
		static auto A1_high_y_stat_hist_ = (*file.Get<TH1F>("Table 2/Hist1D_y2_e1"));
		static auto A2_high_y_stat_hist_ = (*file.Get<TH1F>("Table 2/Hist1D_y3_e1"));

		// syst
		static auto A0_low_y_syst_hist_ = (*file.Get<TH1F>("Table 1/Hist1D_y1_e2"));
		static auto A1_low_y_syst_hist_ = (*file.Get<TH1F>("Table 1/Hist1D_y2_e2"));
		static auto A2_low_y_syst_hist_ = (*file.Get<TH1F>("Table 1/Hist1D_y3_e2"));
		static auto A0_high_y_syst_hist_ = (*file.Get<TH1F>("Table 2/Hist1D_y1_e2"));
		static auto A1_high_y_syst_hist_ = (*file.Get<TH1F>("Table 2/Hist1D_y2_e2"));
		static auto A2_high_y_syst_hist_ = (*file.Get<TH1F>("Table 2/Hist1D_y3_e2"));

		static auto coefficients_hists = CoefficientsHists();
		coefficients_hists.A0_low_y = static_cast<void *>(&A0_low_y_hist_);
		coefficients_hists.A1_low_y = static_cast<void *>(&A1_low_y_hist_);
		coefficients_hists.A2_low_y = static_cast<void *>(&A2_low_y_hist_);
		coefficients_hists.A0_high_y = static_cast<void *>(&A0_high_y_hist_);
		coefficients_hists.A1_high_y = static_cast<void *>(&A1_high_y_hist_);
		coefficients_hists.A2_high_y = static_cast<void *>(&A2_high_y_hist_);
		coefficients_hists.A0_low_y_stat = static_cast<void *>(&A0_low_y_stat_hist_);
		coefficients_hists.A1_low_y_stat = static_cast<void *>(&A1_low_y_stat_hist_);
		coefficients_hists.A2_low_y_stat = static_cast<void *>(&A2_low_y_stat_hist_);
		coefficients_hists.A0_high_y_stat = static_cast<void *>(&A0_high_y_stat_hist_);
		coefficients_hists.A1_high_y_stat = static_cast<void *>(&A1_high_y_stat_hist_);
		coefficients_hists.A2_high_y_stat = static_cast<void *>(&A2_high_y_stat_hist_);
		coefficients_hists.A0_low_y_syst = static_cast<void *>(&A0_low_y_syst_hist_);
		coefficients_hists.A1_low_y_syst = static_cast<void *>(&A1_low_y_syst_hist_);
		coefficients_hists.A2_low_y_syst = static_cast<void *>(&A2_low_y_syst_hist_);
		coefficients_hists.A0_high_y_syst = static_cast<void *>(&A0_high_y_syst_hist_);
		coefficients_hists.A1_high_y_syst = static_cast<void *>(&A1_high_y_syst_hist_);
		coefficients_hists.A2_high_y_syst = static_cast<void *>(&A2_high_y_syst_hist_);

		return static_cast<void *>(&coefficients_hists);
	}

	typedef struct AngularCoefficients
	{
		double A0;
		double A1;
		double A2;
	} AngularCoefficients;

	AngularCoefficients get_angular_coefficients(const CoefficientsHists &coefficients_hists, const double &z_rapidity, const double &qT, const int &syst)
	{
		auto angular_coefficients = AngularCoefficients();
		if (fabs(z_rapidity) < 1.0)
		{
			const auto h_A0 = static_cast<TH1F *>(coefficients_hists.A0_low_y);
			const auto h_A1 = static_cast<TH1F *>(coefficients_hists.A1_low_y);
			const auto h_A2 = static_cast<TH1F *>(coefficients_hists.A2_low_y);
			angular_coefficients.A0 = get_bin_content(qT, h_A0);
			angular_coefficients.A1 = get_bin_content(qT, h_A1);
			angular_coefficients.A2 = get_bin_content(qT, h_A2);
			if (syst == -1 || syst == +1)
			{
				const auto h_A0_stat = static_cast<TH1F *>(coefficients_hists.A0_low_y_stat);
				const auto h_A1_stat = static_cast<TH1F *>(coefficients_hists.A1_low_y_stat);
				const auto h_A2_stat = static_cast<TH1F *>(coefficients_hists.A2_low_y_stat);
				const auto h_A0_syst = static_cast<TH1F *>(coefficients_hists.A0_low_y_syst);
				const auto h_A1_syst = static_cast<TH1F *>(coefficients_hists.A1_low_y_syst);
				const auto h_A2_syst = static_cast<TH1F *>(coefficients_hists.A2_low_y_syst);
				angular_coefficients.A0 += ((double)syst) * sqrt(pow2(get_bin_content(qT, h_A0_stat)) + pow2(get_bin_content(qT, h_A0_syst)));
				angular_coefficients.A1 += ((double)syst) * sqrt(pow2(get_bin_content(qT, h_A1_stat)) + pow2(get_bin_content(qT, h_A1_syst)));
				angular_coefficients.A2 += ((double)syst) * sqrt(pow2(get_bin_content(qT, h_A2_stat)) + pow2(get_bin_content(qT, h_A2_syst)));
			}
		}
		else
		{
			const auto h_A0 = static_cast<TH1F *>(coefficients_hists.A0_high_y);
			const auto h_A1 = static_cast<TH1F *>(coefficients_hists.A1_high_y);
			const auto h_A2 = static_cast<TH1F *>(coefficients_hists.A2_high_y);
			angular_coefficients.A0 = get_bin_content(qT, h_A0);
			angular_coefficients.A1 = get_bin_content(qT, h_A1);
			angular_coefficients.A2 = get_bin_content(qT, h_A2);
			if (syst == -1 || syst == +1)
			{
				const auto h_A0_stat = static_cast<TH1F *>(coefficients_hists.A0_high_y_stat);
				const auto h_A1_stat = static_cast<TH1F *>(coefficients_hists.A1_high_y_stat);
				const auto h_A2_stat = static_cast<TH1F *>(coefficients_hists.A2_high_y_stat);
				const auto h_A0_syst = static_cast<TH1F *>(coefficients_hists.A0_high_y_syst);
				const auto h_A1_syst = static_cast<TH1F *>(coefficients_hists.A1_high_y_syst);
				const auto h_A2_syst = static_cast<TH1F *>(coefficients_hists.A2_high_y_syst);
				angular_coefficients.A0 += ((double)syst) * sqrt(pow2(get_bin_content(qT, h_A0_stat)) + pow2(get_bin_content(qT, h_A0_syst)));
				angular_coefficients.A1 += ((double)syst) * sqrt(pow2(get_bin_content(qT, h_A1_stat)) + pow2(get_bin_content(qT, h_A1_syst)));
				angular_coefficients.A2 += ((double)syst) * sqrt(pow2(get_bin_content(qT, h_A2_stat)) + pow2(get_bin_content(qT, h_A2_syst)));
			}
		}

		return angular_coefficients;
	}

	double get_weight_full_method(double pt_mu_plus,
								  double eta_mu_plus,
								  double phi_mu_plus,
								  double pt_mu_minus,
								  double eta_mu_minus,
								  double phi_mu_minus,
								  double pt_photon,
								  double eta_photon,
								  double phi_photon,
								  void *coefficients_hists_ptr,
								  double z_rapidity,
								  double qT,
								  int syst)
	{
		// check syst input
		if (syst != -1 && syst != 0 && syst != +1)
		{
			printf("[ERROR] [POLARIZATION WEIGHT] Invalid argument. \"syst\" should be -1 (DOWN), 0 (NOMINAL) or +1 (UP).");
			exit(-1);
		}

		const auto coefficients_hists = *(static_cast<CoefficientsHists *>(coefficients_hists_ptr));

		const auto cs_angles = get_CSAngles(pt_mu_plus, eta_mu_plus, phi_mu_plus, pt_mu_minus, eta_mu_minus, phi_mu_minus, pt_photon, eta_photon, phi_photon);
		const double T = cs_angles.BIG_THETA;
		const double P = cs_angles.BIG_PHI;
		const double t = cs_angles.theta;
		const double p = cs_angles.phi;

		const auto angular_coefficients = get_angular_coefficients(coefficients_hists, z_rapidity, qT, syst);
		const double A0 = angular_coefficients.A0;
		const double A1 = angular_coefficients.A1;
		const double A2 = angular_coefficients.A2;

		double w = 0.0;
		w += -4;
		w += A0 * (1 - cos2(T) - cos2(t));
		w += (4 - 3 * A0) * cos2(T) * cos2(t);
		w += A0 * sin2(T) * sin2(t) * cos(2 * (P - p));
		w += 0.5 * (2 - A0) * sin(2 * T) * sin(2 * t) * cos(P - p);
		w += 0.5 * A2 * (sin(2 * T) * sin(2 * t) * cos(P + p) + 2 * sin2(T) * sin2(t) * (cos(2 * P) + cos(2 * p)));
		w += A1 * ((1 + cos2(T)) * sin2(t) * cos(p) + (1 + cos2(t)) * sin(2 * T) * cos(P) + sin2(T) * sin(2 * t) * cos(2 * P - p) + sin2(t) * sin(2 * T) * cos(2 * p - P));

		return w;
	}
}
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
	const double MUON_MASS = 0.1056583755;
	const double PROTON_MASS = 0.93827208816;
	const double BEAM_ENERGY = 13000.0;

	double upsilon_polarization_cos_angle(double pt_mu_plus, double eta_mu_plus, double phi_mu_plus, double pt_mu_minus, double eta_mu_minus, double phi_mu_minus)
	{
		auto gen_mu_plus = PtEtaPhiMVector(pt_mu_plus, eta_mu_plus, phi_mu_plus, MUON_MASS);
		auto gen_mu_minus = PtEtaPhiMVector(pt_mu_minus, eta_mu_minus, phi_mu_minus, MUON_MASS);

		auto upsilon_CM = (gen_mu_plus + gen_mu_minus);

		auto gen_mu_plus_upsilon_frame = VectorUtil::boost(gen_mu_plus, upsilon_CM.BoostToCM());

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
		auto polWgt_Nominal = 1.0;								  // unpolarized
		auto polWgt_PLUS = (3. / 4.) * (1. + pow(cos_angle, 2));  // transverse
		auto polWgt_MINUS = (3. / 2.) * (1. - pow(cos_angle, 2)); // longitudinal

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
		auto gen_mu_plus = PtEtaPhiMVector(pt_mu_plus, eta_mu_plus, phi_mu_plus, MUON_MASS);
		auto gen_mu_minus = PtEtaPhiMVector(pt_mu_minus, eta_mu_minus, phi_mu_minus, MUON_MASS);
		auto upsilon = (gen_mu_plus + gen_mu_minus);
		auto photon = PtEtaPhiMVector(pt_photon, eta_photon, phi_photon, 0.0);
		auto boson = (upsilon + photon);

		// define beam projectile and target
		const double BEAM_MOMEMTUM = sqrt(pow2(BEAM_ENERGY) - pow2(PROTON_MASS));

		PxPyPzEVector projectile; // projectile
		PxPyPzEVector target;	  // target

		if (boson.pz() >= 0)
		{
			projectile = PxPyPzEVector(0., 0., BEAM_MOMEMTUM, BEAM_ENERGY); // projectile
			target = PxPyPzEVector(0., 0., -BEAM_MOMEMTUM, BEAM_ENERGY);	// target
		}
		else
		{
			projectile = PxPyPzEVector(0., 0., -BEAM_MOMEMTUM, BEAM_ENERGY); // projectile
			target = PxPyPzEVector(0., 0., BEAM_MOMEMTUM, BEAM_ENERGY);		 // target
		}

		// beta to boost to CS frame
		auto beta_upsilon_photon = (upsilon + photon).BoostToCM();

		// boost to CS frame
		auto upsilon_CS_frame = VectorUtil::boost(upsilon, beta_upsilon_photon);
		// auto photon_CS_frame = VectorUtil::boost(photon, beta_upsilon_photon); // not used
		auto projectile_CS_frame = VectorUtil::boost(projectile, beta_upsilon_photon);
		auto target_CS_frame = VectorUtil::boost(target, beta_upsilon_photon);

		// Determine x, y and z axis for the CS frame
		auto z_axis_CS_frame = (projectile_CS_frame.Vect().Unit() - target_CS_frame.Vect().Unit()).Unit();
		auto y_axis_CS_frame = (projectile_CS_frame.Vect().Unit().Cross(target_CS_frame.Vect().Unit())).Unit();
		auto x_axis_CS_frame = y_axis_CS_frame.Cross(z_axis_CS_frame).Unit();

		// BIG_PHI
		double BIG_PHI = TMath::ATan2((upsilon_CS_frame.Vect()).Dot(y_axis_CS_frame), ((upsilon_CS_frame.Vect()).Dot(x_axis_CS_frame)));

		// cos_BIG_THETA
		double cos_BIG_THETA = z_axis_CS_frame.Dot((upsilon_CS_frame.Vect()).Unit());
		// Theta CS is not properly defined for Like-Sign muons
		if (cos_BIG_THETA < 0)
			cos_BIG_THETA = -cos_BIG_THETA;
		double BIG_THETA = acos(cos_BIG_THETA);

		// rigid translation to the upsilon CM frame
		// z axis is preserved
		auto gen_mu_plus_upsilon_frame = VectorUtil::boost(gen_mu_plus, upsilon.BoostToCM());
		// auto gen_mu_minus_upsilon_frame = VectorUtil::boost(gen_mu_minus, upsilon.BoostToCM()); // not used

		// phi
		double phi = TMath::ATan2((gen_mu_plus_upsilon_frame.Vect()).Dot(y_axis_CS_frame), ((gen_mu_plus_upsilon_frame.Vect()).Dot(x_axis_CS_frame)));

		// cos_theta
		double cos_theta = z_axis_CS_frame.Dot((gen_mu_plus_upsilon_frame.Vect()).Unit());
		// Theta CS is not properly defined for Like-Sign muons
		if (cos_theta < 0)
			cos_theta = -cos_theta;
		double theta = acos(cos_theta);

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

	// void free_CoefficientsHists(CoefficientsHists *p)
	// {
	// 	free(p);
	// }

	CoefficientsHists load_CMS_data()
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

		auto coefficients_hists = CoefficientsHists();
		coefficients_hists.A0_low_y = static_cast<void*>(&A0_low_y_hist_);
		coefficients_hists.A1_low_y = static_cast<void*>(&A1_low_y_hist_);
		coefficients_hists.A2_low_y = static_cast<void*>(&A2_low_y_hist_);
		coefficients_hists.A0_high_y = static_cast<void*>(&A0_high_y_hist_);
		coefficients_hists.A1_high_y = static_cast<void*>(&A1_high_y_hist_);
		coefficients_hists.A2_high_y = static_cast<void*>(&A2_high_y_hist_);
		coefficients_hists.A0_low_y_stat = static_cast<void*>(&A0_low_y_stat_hist_);
		coefficients_hists.A1_low_y_stat = static_cast<void*>(&A1_low_y_stat_hist_);
		coefficients_hists.A2_low_y_stat = static_cast<void*>(&A2_low_y_stat_hist_);
		coefficients_hists.A0_high_y_stat = static_cast<void*>(&A0_high_y_stat_hist_);
		coefficients_hists.A1_high_y_stat = static_cast<void*>(&A1_high_y_stat_hist_);
		coefficients_hists.A2_high_y_stat = static_cast<void*>(&A2_high_y_stat_hist_);
		coefficients_hists.A0_low_y_syst = static_cast<void*>(&A0_low_y_syst_hist_);
		coefficients_hists.A1_low_y_syst = static_cast<void*>(&A1_low_y_syst_hist_);
		coefficients_hists.A2_low_y_syst = static_cast<void*>(&A2_low_y_syst_hist_);
		coefficients_hists.A0_high_y_syst = static_cast<void*>(&A0_high_y_syst_hist_);
		coefficients_hists.A1_high_y_syst = static_cast<void*>(&A1_high_y_syst_hist_);
		coefficients_hists.A2_high_y_syst = static_cast<void*>(&A2_high_y_syst_hist_);

		return coefficients_hists;
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
			auto h_A0 = static_cast<TH1F *>(coefficients_hists.A0_low_y);
			auto h_A1 = static_cast<TH1F *>(coefficients_hists.A1_low_y);
			auto h_A2 = static_cast<TH1F *>(coefficients_hists.A2_low_y);
			angular_coefficients.A0 = get_bin_content(qT, h_A0);
			angular_coefficients.A1 = get_bin_content(qT, h_A1);
			angular_coefficients.A2 = get_bin_content(qT, h_A2);
			if (syst == -1 || syst == +1)
			{
				auto h_A0_stat = static_cast<TH1F *>(coefficients_hists.A0_low_y_stat);
				auto h_A1_stat = static_cast<TH1F *>(coefficients_hists.A1_low_y_stat);
				auto h_A2_stat = static_cast<TH1F *>(coefficients_hists.A2_low_y_stat);
				auto h_A0_syst = static_cast<TH1F *>(coefficients_hists.A0_low_y_syst);
				auto h_A1_syst = static_cast<TH1F *>(coefficients_hists.A1_low_y_syst);
				auto h_A2_syst = static_cast<TH1F *>(coefficients_hists.A2_low_y_syst);
				angular_coefficients.A0 += ((double)syst) * sqrt(pow2(get_bin_content(qT, h_A0_stat)) + pow2(get_bin_content(qT, h_A0_syst)));
				angular_coefficients.A1 += ((double)syst) * sqrt(pow2(get_bin_content(qT, h_A1_stat)) + pow2(get_bin_content(qT, h_A1_syst)));
				angular_coefficients.A2 += ((double)syst) * sqrt(pow2(get_bin_content(qT, h_A2_stat)) + pow2(get_bin_content(qT, h_A2_syst)));
			}
		}
		else
		{
			auto h_A0 = static_cast<TH1F *>(coefficients_hists.A0_high_y);
			auto h_A1 = static_cast<TH1F *>(coefficients_hists.A1_high_y);
			auto h_A2 = static_cast<TH1F *>(coefficients_hists.A2_high_y);
			angular_coefficients.A0 = get_bin_content(qT, h_A0);
			angular_coefficients.A1 = get_bin_content(qT, h_A1);
			angular_coefficients.A2 = get_bin_content(qT, h_A2);
			if (syst == -1 || syst == +1)
			{
				auto h_A0_stat = static_cast<TH1F *>(coefficients_hists.A0_high_y_stat);
				auto h_A1_stat = static_cast<TH1F *>(coefficients_hists.A1_high_y_stat);
				auto h_A2_stat = static_cast<TH1F *>(coefficients_hists.A2_high_y_stat);
				auto h_A0_syst = static_cast<TH1F *>(coefficients_hists.A0_high_y_syst);
				auto h_A1_syst = static_cast<TH1F *>(coefficients_hists.A1_high_y_syst);
				auto h_A2_syst = static_cast<TH1F *>(coefficients_hists.A2_high_y_syst);
				angular_coefficients.A0 += ((double)syst) * sqrt(pow2(get_bin_content(qT, h_A0_stat)) + pow2(get_bin_content(qT, h_A0_syst)));
				angular_coefficients.A1 += ((double)syst) * sqrt(pow2(get_bin_content(qT, h_A1_stat)) + pow2(get_bin_content(qT, h_A1_syst)));
				angular_coefficients.A2 += ((double)syst) * sqrt(pow2(get_bin_content(qT, h_A2_stat)) + pow2(get_bin_content(qT, h_A2_syst)));
			}
		}

		return angular_coefficients;
	}

	// How to pass string from ctpes:
	// https://stackoverflow.com/questions/37966432/how-to-pass-const-char-from-python-to-c-function
	double get_weight_full_method(double pt_mu_plus,
								  double eta_mu_plus,
								  double phi_mu_plus,
								  double pt_mu_minus,
								  double eta_mu_minus,
								  double phi_mu_minus,
								  double pt_photon,
								  double eta_photon,
								  double phi_photon,
								  const CoefficientsHists &coefficients_hists,
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

		auto cs_angles = get_CSAngles(pt_mu_plus, eta_mu_plus, phi_mu_plus, pt_mu_minus, eta_mu_minus, phi_mu_minus, pt_photon, eta_photon, phi_photon);
		double T = cs_angles.BIG_THETA;
		double P = cs_angles.BIG_PHI;
		double t = cs_angles.theta;
		double p = cs_angles.phi;

		auto angular_coefficients = get_angular_coefficients(coefficients_hists, z_rapidity, qT, syst);
		double A0 = angular_coefficients.A0;
		double A1 = angular_coefficients.A1;
		double A2 = angular_coefficients.A2;

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
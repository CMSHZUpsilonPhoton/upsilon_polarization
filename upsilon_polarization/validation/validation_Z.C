
#include <any>
#include <map>
#include <string>

#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "ROOT/RVec.hxx"
#include "TDatabasePDG.h"
#include "TH1F.h"
#include "TMath.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "Math/VectorUtil.h"

#include "../src/upsilon_polarization.cxx"

#include "tools/tabulate.hpp"
using namespace tabulate;

using ROOT::VecOps::RVec;

void fill_extreme_scenarios(const ROOT::Math::PtEtaPhiMVector &gamma,
                            const ROOT::Math::PtEtaPhiMVector &mu_plus,
                            const ROOT::Math::PtEtaPhiMVector &mu_minus,
                            const ROOT::Math::PtEtaPhiMVector &upsilon,
                            const ROOT::Math::PtEtaPhiMVector &boson,
                            const ROOT::Math::PxPyPzEVector &proton_plus,
                            const ROOT::Math::PxPyPzEVector &proton_minus,
                            TH1F &h_cos_BIG_THETA_unweighted,
                            TH1F &h_BIG_PHI_unweighted,
                            TH1F &h_cos_theta_unweighted,
                            TH1F &h_phi_unweighted,
                            TH1F &h_cos_BIG_THETA_weighted,
                            TH1F &h_BIG_PHI_weighted,
                            TH1F &h_cos_theta_weighted,
                            TH1F &h_phi_weighted)
{
    PxPyPzEVector beam;   // beam
    PxPyPzEVector target; // target

    // get beam x target
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

    // HZ rest frame
    const auto upsilon_HZ_rest_frame = VectorUtil::boost(upsilon, boson.BoostToCM());
    const auto gamma_HZ_rest_frame = VectorUtil::boost(gamma, boson.BoostToCM());
    const auto mu_plus_HZ_rest_frame = VectorUtil::boost(mu_plus, boson.BoostToCM());
    const auto mu_minus_HZ_rest_frame = VectorUtil::boost(mu_minus, boson.BoostToCM());
    const auto beam_HZ_rest_frame = VectorUtil::boost(beam, boson.BoostToCM());

    auto cos_BIG_THETA = upsilon_HZ_rest_frame.Vect().Unit().Dot(beam_HZ_rest_frame.Vect().Unit());

    auto n_beam_gamma_plane = beam_HZ_rest_frame.Vect().Unit().Cross(gamma_HZ_rest_frame.Vect().Unit());
    auto n_dimuon_plane = mu_plus_HZ_rest_frame.Vect().Unit().Cross(mu_minus_HZ_rest_frame.Vect().Unit());
    auto BIG_PHI = acos(n_beam_gamma_plane.Dot(n_dimuon_plane));

    // Upsilon rest frame
    const auto upsilon_Upsilon_rest_frame = VectorUtil::boost(upsilon, upsilon.BoostToCM());
    const auto gamma_Upsilon_rest_frame = VectorUtil::boost(gamma, upsilon.BoostToCM());
    const auto mu_plus_Upsilon_rest_frame = VectorUtil::boost(mu_plus, upsilon.BoostToCM());
    const auto mu_minus_Upsilon_rest_frame = VectorUtil::boost(mu_minus, upsilon.BoostToCM());
    const auto beam_Upsilon_rest_frame = VectorUtil::boost(beam, upsilon.BoostToCM());

    auto cos_theta = mu_plus_Upsilon_rest_frame.Vect().Unit().Dot(upsilon.Vect().Unit());
    auto phi = 1.0; // does not apply in this scenario (?)

    h_cos_BIG_THETA_unweighted.Fill(cos_BIG_THETA);
    h_BIG_PHI_unweighted.Fill(BIG_PHI);
    h_cos_theta_unweighted.Fill(cos_theta);
    h_phi_unweighted.Fill(phi);
    h_cos_BIG_THETA_weighted.Fill(cos_BIG_THETA);
    h_BIG_PHI_weighted.Fill(BIG_PHI);
    h_cos_theta_weighted.Fill(cos_theta);
    h_phi_weighted.Fill(phi);
}
void fill_alternative() {}
void fill_full() {}

bool get_bit(int b, int bitNumber)
{
    return (b & (1 << bitNumber)) != 0;
}

// HEP Object
class HEPObjectVec
{
private:
    const size_t size_;
    std::map<std::string, std::any> features;

public:
    HEPObjectVec(std::map<std::string, std::any> _features, size_t _size = 1) : size_(_size), features(_features)
    {
        if (size_ < 1)
        {
            throw std::invalid_argument("Size should be > 0.");
        }
    }

    void add(const std::string name, std::any feature)
    {
        features[name] = feature;
    }

    template <typename T>
    T val(const std::string feature)
    {
        return std::any_cast<T>(features[feature]);
    }

    template <typename T>
    T value(const std::string feature)
    {
        return val<T>(feature);
    }

    template <typename T>
    RVec<T> vec(const std::string feature)
    {
        return val<RVec<T>>(feature);
    }

    template <typename T>
    RVec<T> vector(const std::string feature)
    {
        return vec<T>(feature);
    }

    size_t size()
    {
        return size_;
    }
};

auto db = TDatabasePDG();

void validation_Z(int method)
{
    // gSystem->Load("upsilon_polarization/lib/upsilon_polarization");
    // std::cout << get_weight_extreme_scenarios(0.5, -1) << std::endl;

    ROOT::RDataFrame nanoaod_df("Events", "validation_data/ZToUpsilon1SGamma_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL18NanoAODv9.root");

    // nanoaod_df.Describe().Print();

    // define histograms
    auto n_bins_cos = 20;
    auto n_bins_phi = 20;

    auto h_cos_BIG_THETA_unweighted = TH1F("h_cos_BIG_THETA_unweighted", "cos_BIG_THETA", n_bins_cos, -1.0, 1.0);
    h_cos_BIG_THETA_unweighted.Sumw2();

    auto h_BIG_PHI_unweighted = TH1F("h_BIG_PHI_unweighted", "BIG_PHI", n_bins_phi, -TMath::Pi(), TMath::Pi());
    h_BIG_PHI_unweighted.Sumw2();

    auto h_cos_theta_unweighted = TH1F("h_cos_theta_unweighted", "cos_theta", n_bins_cos, -1.0, 1.0);
    h_cos_theta_unweighted.Sumw2();

    auto h_phi_unweighted = TH1F("h_phi_unweighted", "phi", n_bins_phi, -TMath::Pi(), TMath::Pi());
    h_phi_unweighted.Sumw2();

    auto h_cos_BIG_THETA_weighted = TH1F("h_cos_BIG_THETA_weighted", "cos_BIG_THETA", n_bins_cos, -1.0, 1.0);
    h_cos_BIG_THETA_weighted.Sumw2();

    auto h_BIG_PHI_weighted = TH1F("h_BIG_PHI_weighted", "BIG_PHI", n_bins_phi, -TMath::Pi(), TMath::Pi());
    h_BIG_PHI_weighted.Sumw2();

    auto h_cos_theta_weighted = TH1F("h_cos_theta_weighted", "cos_theta", n_bins_cos, -1.0, 1.0);
    h_cos_theta_weighted.Sumw2();

    auto h_phi_weighted = TH1F("h_phi_weighted", "phi", n_bins_phi, -TMath::Pi(), TMath::Pi());
    h_phi_weighted.Sumw2();

    auto event_idx = 0;
    std::cout << "\n\nStarting event loop ..." << std::endl;
    nanoaod_df.Foreach(
        [&](UInt_t nGenPart,
            RVec<Float_t> GenPart_eta,
            RVec<Int_t> GenPart_genPartIdxMother,
            RVec<Float_t> GenPart_mass,
            RVec<Int_t> GenPart_pdgId,
            RVec<Float_t> GenPart_phi,
            RVec<Float_t> GenPart_pt,
            RVec<Int_t> GenPart_status,
            RVec<Int_t> GenPart_statusFlags)
        {
            auto gen_part = HEPObjectVec({
                                             {"eta", std::any(GenPart_eta)},
                                             {"genPartIdxMother", std::any(GenPart_genPartIdxMother)},
                                             {"mass", std::any(GenPart_mass)},
                                             {"pdgId", std::any(GenPart_pdgId)},
                                             {"phi", std::any(GenPart_phi)},
                                             {"pt", std::any(GenPart_pt)},
                                             {"status", std::any(GenPart_status)},
                                             {"statusFlags", std::any(GenPart_statusFlags)},
                                         },
                                         nGenPart);

            event_idx++;

            // loop over particles
            auto idx_z0 = -1;
            auto idx_upsilon = -1;
            auto idx_gamma = -1;
            auto idx_mu_plus = -1;
            auto idx_mu_minus = -1;
            for (unsigned int idx_part = 0; idx_part < gen_part.size(); idx_part++)
            {
                if (gen_part.vector<int>("pdgId")[idx_part] == 23)
                {
                    idx_z0 = idx_part;
                }
                if (gen_part.vector<int>("pdgId")[idx_part] == 553 && (gen_part.vector<int>("genPartIdxMother")[idx_part] == idx_z0 || gen_part.vector<int>("genPartIdxMother")[idx_part] == idx_upsilon))
                {
                    idx_upsilon = idx_part;
                }
                if (gen_part.vector<int>("pdgId")[idx_part] == 22 && (gen_part.vector<int>("genPartIdxMother")[idx_part] == idx_z0 || gen_part.vector<int>("genPartIdxMother")[idx_part] == idx_gamma))
                {
                    idx_gamma = idx_part;
                }
                if (gen_part.vector<int>("pdgId")[idx_part] == 13 && gen_part.vector<int>("genPartIdxMother")[idx_part] == idx_upsilon)
                {
                    idx_mu_plus = idx_part;
                }
                if (gen_part.vector<int>("pdgId")[idx_part] == -13 && gen_part.vector<int>("genPartIdxMother")[idx_part] == idx_upsilon)
                {
                    idx_mu_minus = idx_part;
                }
            }

            if (ROOT::VecOps::Any(RVec<int>{idx_z0, idx_upsilon, idx_gamma, idx_mu_plus, idx_mu_minus} < 0))
            {
                throw std::runtime_error("Problem looking for particle indexes.\nRVec<int>{idx_z0, idx_upsilon, idx_gamma, idx_mu_plus, idx_mu_minus}");
            }

            // build particle lorentz vectors
            const auto gen_part_pt = gen_part.vector<float>("pt");
            const auto gen_part_eta = gen_part.vector<float>("eta");
            const auto gen_part_phi = gen_part.vector<float>("phi");
            const auto gen_part_mass = gen_part.vector<float>("mass");

            const auto _boson = PtEtaPhiMVector(
                gen_part_pt[idx_z0],
                gen_part_eta[idx_z0],
                gen_part_phi[idx_z0],
                gen_part_mass[idx_z0]);
            const auto _upsilon = PtEtaPhiMVector(
                gen_part_pt[idx_upsilon],
                gen_part_eta[idx_upsilon],
                gen_part_phi[idx_upsilon],
                gen_part_mass[idx_upsilon]);
            const auto gamma = PtEtaPhiMVector(
                gen_part_pt[idx_gamma],
                gen_part_eta[idx_gamma],
                gen_part_phi[idx_gamma],
                0.0);
            const auto mu_plus = PtEtaPhiMVector(
                gen_part_pt[idx_mu_plus],
                gen_part_eta[idx_mu_plus],
                gen_part_phi[idx_mu_plus],
                // 0);
                MUON_MASS);
            const auto mu_minus = PtEtaPhiMVector(
                gen_part_pt[idx_mu_minus],
                gen_part_eta[idx_mu_minus],
                gen_part_phi[idx_mu_minus],
                // 0);
                MUON_MASS);

            const auto upsilon = mu_plus + mu_minus;
            const auto boson = mu_plus + mu_minus + gamma;
            // const auto upsilon = _upsilon;
            // const auto boson = _boson;

            
            // define beam beam and target
            const double PROTON_MOMEMTUM = sqrt(pow2(PROTON_ENERGY) - pow2(PROTON_MASS));
            const auto proton_plus = PxPyPzEVector(0., 0., PROTON_MOMEMTUM, PROTON_ENERGY);
            const auto proton_minus = PxPyPzEVector(0., 0., -PROTON_MOMEMTUM, PROTON_ENERGY);

            if (method == 1)
            {
                fill_extreme_scenarios(gamma, mu_plus, mu_minus, upsilon, boson, proton_plus, proton_minus, h_cos_BIG_THETA_unweighted, h_BIG_PHI_unweighted, h_cos_theta_unweighted, h_phi_unweighted, h_cos_BIG_THETA_weighted, h_BIG_PHI_weighted, h_cos_theta_weighted, h_phi_weighted);
            }
            if (method == 2)
            {
                fill_alternative();
            }
            if (method == 3)
            {
                fill_full();
            }
            if (event_idx % 100000 == 0)
            {
                std::cout << "--> Processed events: " << event_idx << " ..." << std::endl;

                    std::cout << " " << std::endl;
                    std::cout << " " << std::endl;
                    auto particles = Table();
                    particles.add_row({"Index",
                                       "Particle",
                                       "PDG ID",
                                       "Status",
                                       "Status Flags",
                                       "Mother Index",
                                       "Mother Name",
                                       "pt",
                                       "eta",
                                       "phi",
                                       "mass"});
                    // center-align and color header cells
                    for (unsigned int i = 0; i < particles[0].size(); i++)
                    {
                        particles[0][i].format().font_color(Color::yellow).font_align(FontAlign::center).font_style({FontStyle::bold});
                    }

                    for (unsigned int i = 0; i < gen_part.size(); i++)
                    {
                        auto pdgId = gen_part.vector<int>("pdgId")[i];
                        auto status = gen_part.vector<int>("status")[i];
                        auto statusFlags = gen_part.vector<int>("statusFlags")[i];
                        auto IdxMother = gen_part.vector<int>("genPartIdxMother")[i];
                        auto pt = gen_part.vector<float>("pt")[i];
                        auto eta = gen_part.vector<float>("eta")[i];
                        auto phi = gen_part.vector<float>("phi")[i];
                        auto mass = gen_part.vector<float>("mass")[i];
                        particles.add_row({std::to_string(i),
                                           db.GetParticle(pdgId)->GetName(),
                                           std::to_string(pdgId),
                                           std::to_string(status),
                                           std::to_string(statusFlags),
                                           std::to_string(IdxMother),
                                           db.GetParticle(GenPart_pdgId[IdxMother])->GetName(),
                                           std::to_string(pt),
                                           std::to_string(eta),
                                           std::to_string(phi),
                                           std::to_string(mass)});
                    }
                    std::cout << particles << std::endl;
                    std::cout << " " << std::endl;
                    std::cout << " " << std::endl;
                std::cout << "Print properties:" << std::endl;
                std::cout << "Boson pt: " << _boson.pt() << std::endl;
                std::cout << "Boson eta: " << _boson.eta() << std::endl;
                std::cout << "Boson phi: " << _boson.phi() << std::endl;
                std::cout << "Boson mass: " << _boson.mass() << std::endl;
                std::cout << "Upsilon pt: " << _upsilon.pt() << std::endl;
                std::cout << "Upsilon eta: " << _upsilon.eta() << std::endl;
                std::cout << "Upsilon phi: " << _upsilon.phi() << std::endl;
                std::cout << "Upsilon mass: " << _upsilon.mass() << std::endl;
                std::cout << "************************************************" << std::endl;

                std::cout << "Print diff:" << std::endl;
                std::cout << "Boson pt: " << _boson.pt() - boson.pt() << std::endl;
                std::cout << "Boson eta: " << _boson.eta() - boson.eta() << std::endl;
                std::cout << "Boson phi: " << _boson.phi() - boson.phi() << std::endl;
                std::cout << "Boson mass: " << _boson.mass() - boson.mass() << std::endl;
                std::cout << "Upsilon pt: " << _upsilon.pt() - upsilon.pt() << std::endl;
                std::cout << "Upsilon eta: " << _upsilon.eta() - upsilon.eta() << std::endl;
                std::cout << "Upsilon phi: " << _upsilon.phi() - upsilon.phi() << std::endl;
                std::cout << "Upsilon mass: " << _upsilon.mass() - upsilon.mass() << std::endl;
                std::cout << "************************************************" << std::endl;
            }
        },
        // columns names
        {"nGenPart",
         "GenPart_eta",
         "GenPart_genPartIdxMother",
         "GenPart_mass",
         "GenPart_pdgId",
         "GenPart_phi",
         "GenPart_pt",
         "GenPart_status",
         "GenPart_statusFlags"});

    // normalize histograms to unit
    h_cos_BIG_THETA_unweighted.Scale(1.0 / h_cos_BIG_THETA_unweighted.Integral());
    h_BIG_PHI_unweighted.Scale(1.0 / h_BIG_PHI_unweighted.Integral());
    h_cos_theta_unweighted.Scale(1.0 / h_cos_theta_unweighted.Integral());
    h_phi_unweighted.Scale(1.0 / h_phi_unweighted.Integral());
    h_cos_BIG_THETA_unweighted.Scale(1.0 / h_cos_BIG_THETA_unweighted.Integral());
    h_BIG_PHI_unweighted.Scale(1.0 / h_BIG_PHI_unweighted.Integral());
    h_cos_theta_unweighted.Scale(1.0 / h_cos_theta_unweighted.Integral());
    h_phi_unweighted.Scale(1.0 / h_phi_unweighted.Integral());

    // save histograms to file
    std::unique_ptr<TFile> output_file(TFile::Open("validation_outputs/validation_Z.root", "RECREATE"));
    output_file->WriteObject(&h_cos_BIG_THETA_unweighted, "h_cos_BIG_THETA_unweighted");
    output_file->WriteObject(&h_BIG_PHI_unweighted, "h_BIG_PHI_unweighted");
    output_file->WriteObject(&h_cos_theta_unweighted, "h_cos_theta_unweighted");
    output_file->WriteObject(&h_phi_unweighted, "h_phi_unweighted");
    output_file->WriteObject(&h_cos_BIG_THETA_weighted, "h_cos_BIG_THETA_weighted");
    output_file->WriteObject(&h_BIG_PHI_weighted, "h_BIG_PHI_weighted");
    output_file->WriteObject(&h_cos_theta_weighted, "h_cos_theta_weighted");
    output_file->WriteObject(&h_phi_weighted, "h_phi_weighted");

    // h_cos_BIG_THETA_unweighted.Print("all");
    // h_BIG_PHI_unweighted.Print("all");
    // h_cos_theta_unweighted.Print("all");
    // h_phi_unweighted.Print("all");
    // h_cos_BIG_THETA_unweighted.Print("all");
    // h_BIG_PHI_unweighted.Print("all");
    // h_cos_theta_unweighted.Print("all");
    // h_phi_unweighted.Print("all");
}

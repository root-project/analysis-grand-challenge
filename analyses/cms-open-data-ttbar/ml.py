import ROOT
import os
from dataclasses import dataclass

# histogram bin lower limit to use for each ML input feature
bin_low = [0, 0, 0, 0, 50, 50, 50, 50, 25, 25, 25, 25, 0, 0, 0, 0, -1, -1, -1, -1]

# histogram bin upper limit to use for each ML input feature
bin_high = [6, 6, 6, 6, 300, 300, 550, 550, 300, 300, 300, 300, 1, 1, 1, 1, 1, 1, 1, 1]

# names of each ML input feature (used when creating histograms)
feature_names = [
    "deltar_leptonbtoplep", "deltar_w1w2", "deltar_w1btophad", "deltar_w2btophad",
    "mass_leptonbtoplep",   "mass_w1w2",   "mass_w1w2btophad", "pt_w1w2btophad",
    "pt_w1",                "pt_w2",       "pt_btophad",       "pt_btoplep", 
    "btag_w1",              "btag_w2",     "btag_btophad",     "btag_btoplep", 
    "qgl_w1",               "qgl_w2",      "qgl_btophad",      "qgl_btoplep",
]

# labels for each ML input feature (used for plotting)
feature_labels = [
    "Delta R between b_{top-lep} Jet and Lepton",
    "Delta R between the two W Jets",
    "Delta R between first W Jet and b_{top-had} Jet",
    "Delta R between second W Jet and b_{top-had} Jet",
    "Combined Mass of b_{top-lep} Jet and Lepton [GeV]",
    "Combined Mass of the two W Jets [GeV]",
    "Combined Mass of b_{top-had} Jet and the two W Jets [GeV]",
    "Combined p_T of b_{top-had} Jet and the two W Jets [GeV]",
    "p_T of the first W Jet [GeV]",
    "p_T of the second W Jet [GeV]",
    "p_T of the b_{top-had} Jet [GeV]",
    "p_T of the b_{top-lep} Jet [GeV]",
    "btagCSVV2 of the first W Jet",
    "btagCSVV2 of the second W Jet",
    "btagCSVV2 of the b_{top-had} Jet",
    "btagCSVV2 of the b_{top-lep} Jet",
    "Quark vs Gluon likelihood discriminator of the first W Jet",
    "Quark vs Gluon likelihood discriminator of the second W Jet",
    "Quark vs Gluon likelihood discriminator of the b_{top-had} Jet",
    "Quark vs Gluon likelihood discriminator of the b_{top-lep} Jet",
]

@dataclass
class MLHistoConf:
  name: str
  title: str
  binning: (float, float, int) # nbins, low, high

ml_features_config: list[MLHistoConf] = [
    MLHistoConf(name = feature_names[i], title = feature_labels[i], binning = (25, bin_low[i], bin_high[i])) for i in range(len(feature_names))
]




def compile_mlhelpers_cpp(fastforest_path, max_n_jets=6):
    include = os.path.join(fastforest_path, "include")
    lib = os.path.join(fastforest_path, "lib")
    ROOT.gSystem.AddIncludePath(f"-I{include}")
    ROOT.gSystem.AddLinkedLibs(f"-L{lib} -lfastforest")
    ROOT.gSystem.Load(f"{lib}/libfastforest.so.1")
    ROOT.gSystem.CompileMacro("ml_helpers.cpp", "kO")

    ROOT.gInterpreter.Declare(
        """
                              auto models = get_fastforests("models/", 20);
                              auto feven = models["even"];
                              auto fodd = models["odd"];
                              """
    )

    ROOT.gInterpreter.Declare(
        f"""
        size_t max_n_jets = {max_n_jets};
        std::map<int, std::vector<ROOT::RVecI>> permutations = get_permutations_dict(max_n_jets);
        """
    )


def define_lepton_fields(df: ROOT.RDataFrame) -> ROOT.RDataFrame:
    # prepare lepton fields (p4, eta, phi)
    df = (
        df.Define("Electron_phi_cut", "Electron_phi[Electron_mask]")
        .Define("Electron_eta_cut", "Electron_eta[Electron_mask]")
        .Define(
            "Electron_p4",
            """
            ConstructP4 (
                    Electron_pt[Electron_mask], 
                    Electron_eta_cut, 
                    Electron_phi_cut, 
                    Electron_mass[Electron_mask]
            )
            """,
        )
        .Define("Muon_phi_cut", "Muon_phi[Muon_mask]")
        .Define("Muon_eta_cut", "Muon_eta[Muon_mask]")
        .Define(
            "Muon_p4",
            """
            ConstructP4 (
                    Muon_pt[Muon_mask], 
                    Muon_eta_cut, 
                    Muon_phi_cut, 
                    Muon_mass[Muon_mask]
            )            
            """,
        )
        .Define("Lepton_phi", "Concatenate(Electron_phi_cut, Muon_phi_cut)")
        .Define("Lepton_eta", "Concatenate(Electron_eta_cut, Muon_eta_cut)")
        .Define("Lepton_p4", "Concatenate(Electron_p4, Muon_p4)")
    )

    return df


def define_features(df: ROOT.RDataFrame) -> ROOT.RDataFrame:
    # prepare lepton fields (p4, eta, phi)
    df = define_lepton_fields(df)

    # get indexes of four jets
    df = (
        df.Define("W1_idx", "permutations[std::min(Jet_pt_cut.size(), max_n_jets)][0]")
        .Define("W2_idx", "permutations[std::min(Jet_pt_cut.size(), max_n_jets)][1]")
        .Define("bH_idx", "permutations[std::min(Jet_pt_cut.size(), max_n_jets)][2]")
        .Define("bL_idx", "permutations[std::min(Jet_pt_cut.size(), max_n_jets)][3]")
    )

    # # Apply indexes to jets. Jets pt and btagCSVV2 and qgl are features itself (12 features)
    df = (
        df
        # not features themself, but needed to construct features
        .Define("JetW1_phi", "Take(Jet_phi_cut, W1_idx)")
        .Define("JetW2_phi", "Take(Jet_phi_cut, W2_idx)")
        .Define("JetbL_phi", "Take(Jet_phi_cut, bL_idx)")
        .Define("JetbH_phi", "Take(Jet_phi_cut, bH_idx)")
        .Define("JetW1_eta", "Take(Jet_eta_cut, W1_idx)")
        .Define("JetW2_eta", "Take(Jet_eta_cut, W2_idx)")
        .Define("JetbL_eta", "Take(Jet_eta_cut, bL_idx)")
        .Define("JetbH_eta", "Take(Jet_eta_cut, bH_idx)")
        #         # 12 features
        .Define(f"{feature_names[8]}", "Take(Jet_pt_cut, W1_idx)")
        .Define(f"{feature_names[9]}", "Take(Jet_pt_cut, W2_idx)")
        .Define(f"{feature_names[10]}", "Take(Jet_pt_cut, bH_idx)")
        .Define(f"{feature_names[11]}", "Take(Jet_pt_cut, bL_idx)")
        .Define(f"{feature_names[12]}", "Take(Jet_btagCSVV2_cut, W1_idx)")
        .Define(f"{feature_names[13]}", "Take(Jet_btagCSVV2_cut, W2_idx)")
        .Define(f"{feature_names[14]}", "Take(Jet_btagCSVV2_cut, bH_idx)")
        .Define(f"{feature_names[15]}", "Take(Jet_btagCSVV2_cut, bL_idx)")
        .Define("Jet_qgl_cut", "Jet_qgl[Jet_mask]")
        .Define(f"{feature_names[16]}", "Take(Jet_qgl_cut, W1_idx)")
        .Define(f"{feature_names[17]}", "Take(Jet_qgl_cut, W2_idx)")
        .Define(f"{feature_names[18]}", "Take(Jet_qgl_cut, bH_idx)")
        .Define(f"{feature_names[19]}", "Take(Jet_qgl_cut, bL_idx)")
        #        jets' four-momenta
        .Define(
            "JetW1_p4",
            f"ConstructP4({feature_names[8]}, JetW1_eta, JetW1_phi, Take(Jet_mass_cut, W1_idx))",
        )
        .Define(
            "JetW2_p4",
            f"ConstructP4({feature_names[9]}, JetW2_eta, JetW2_phi, Take(Jet_mass_cut, W2_idx))",
        )
        .Define(
            "JetbL_p4",
            f"ConstructP4({feature_names[11]}, JetbL_eta, JetbL_phi, Take(Jet_mass_cut, bL_idx))",
        )
        .Define(
            "JetbH_p4",
            f"ConstructP4({feature_names[10]}, JetbH_eta, JetbH_phi, Take(Jet_mass_cut, bH_idx))",
        )
    )

    # # build features 8 other features
    df = (
        df.Define(
            f"{feature_names[0]}",
            "sqrt(pow(Lepton_eta.at(0)-JetbL_eta,2.)+pow(Lepton_phi.at(0)-JetbL_phi,2.))",
        )
        .Define(
            f"{feature_names[1]}", "sqrt(pow(JetW1_eta-JetW2_eta,2.)+pow(JetW1_phi-JetW2_phi,2.))"
        )
        .Define(
            f"{feature_names[2]}", "sqrt(pow(JetW1_eta-JetbH_eta,2.)+pow(JetW1_phi-JetbH_phi,2.))"
        )
        .Define(
            f"{feature_names[3]}", "sqrt(pow(JetW2_eta-JetbH_eta,2.)+pow(JetW2_phi-JetbH_phi,2.))"
        )
        .Define(
            f"{feature_names[4]}",
            "return Map(JetbL_p4+Lepton_p4.at(0), [] (const ROOT::Math::PxPyPzMVector &p) {return p.M();})",
        )
        .Define(
            f"{feature_names[5]}",
            "return Map(JetW1_p4+JetW2_p4, [] (const ROOT::Math::PxPyPzMVector &p) {return p.M();})",
        )
        .Define(
            f"{feature_names[6]}",
            "return Map(JetW1_p4+JetW2_p4+JetbH_p4, [] (const ROOT::Math::PxPyPzMVector &p) {return p.M();})",
        )
        .Define(
            f"{feature_names[7]}",
            "return Map(JetW1_p4+JetW2_p4+JetbH_p4, [] (const ROOT::Math::PxPyPzMVector &p) {return p.Pt();})",
        )
    )
    # put all features into one collection

    def list_as_str (feature_list: list[str]) -> str:
        """ 
        accept [f0, f1, ..., fN]
        return "{f0, f1, ..., fN}" 
        """
        return feature_list.__str__().replace("[", "{").replace("]", "}").replace("'", "")
    
    df = df.Define("features", f"ROOT::VecOps::RVec<ROOT::RVecF>({list_as_str(feature_names)})")

    return df


def predict_proba(df: ROOT.RDataFrame) -> ROOT.RDataFrame:
    """get probability scores for every permutation in event"""

    # in inference, odd model applied to even events, while even model to odd events
    # read more about inference in dedicated part of AGC documentation: 
    # https://agc.readthedocs.io/en/latest/taskbackground.html#machine-learning-component  
    
    return df.Define(
        "proba",
        """
        bool is_even = (event % 2 == 0);
        const auto& forest = (is_even) ? fodd : feven;
        return inference(features, forest, true);
        """,
    )


def infer_output_ml_features(df: ROOT.RDataFrame) -> ROOT.RDataFrame:
    """
    Choose for each feature the best candidate with the highest probability score.
    Results are features of the best four-jet permutation at each event
    """

    df = predict_proba(df)
    for i, feature in enumerate(feature_names):
        df = df.Define(f"results{i}", f"{feature}[ArgMax(proba)]")

    return df

from dataclasses import dataclass
from typing import Tuple

import ROOT

# histogram bin lower limit to use for each ML input feature
bin_low = [0, 0, 0, 0, 50, 50, 50, 50, 25, 25, 25, 25, 0, 0, 0, 0, -1, -1, -1, -1]

# histogram bin upper limit to use for each ML input feature
bin_high = [6, 6, 6, 6, 300, 300, 550, 550, 300, 300, 300, 300, 1, 1, 1, 1, 1, 1, 1, 1]

# names of each ML input feature (used when creating histograms)
feature_names = [
    "deltar_leptonbtoplep",
    "deltar_w1w2",
    "deltar_w1btophad",
    "deltar_w2btophad",
    "mass_leptonbtoplep",
    "mass_w1w2",
    "mass_w1w2btophad",
    "pt_w1w2btophad",
    "pt_w1",
    "pt_w2",
    "pt_btophad",
    "pt_btoplep",
    "btag_w1",
    "btag_w2",
    "btag_btophad",
    "btag_btoplep",
    "qgl_w1",
    "qgl_w2",
    "qgl_btophad",
    "qgl_btoplep",
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
    binning: Tuple[int, float, float]  # nbins, low, high


ml_features_config: list[MLHistoConf] = [
    MLHistoConf(
        name=feature_names[i], title=feature_labels[i], binning=(25, bin_low[i], bin_high[i])
    )
    for i in range(len(feature_names))
]

def load_cpp(max_n_jets=6):
    # the default value of max_n_jets is the same as in the reference implementation
    # https://github.com/iris-hep/analysis-grand-challenge
    ROOT.gSystem.CompileMacro("ml_helpers.cpp", "kO")

    # Initialize FastForest models.
    # Our BDT models have 20 input features according to the AGC documentation
    # https://agc.readthedocs.io/en/latest/taskbackground.html#machine-learning-component

    ROOT.gInterpreter.Declare(
        # **Conditional directives used to avoid redefinition error during distributed computing**
        # Note:
        # * moving all stuff in `Declare` to `ml_helpers.cpp` cancels the necessity of using `ifndef`
        # * coming soon feature is `gInterpreter.Declare` with automatic header guards
        # https://indico.fnal.gov/event/23628/contributions/240608/attachments/154873/201557/distributed_RDF_padulano_ROOT_workshop_2022.pdf
        """
        #ifndef AGC_MODELS
        #define AGC_MODELS
        
        TMVA::Experimental::RBDT feven("feven", "bdt_even.root");
        TMVA::Experimental::RBDT fodd("fodd", "bdt_odd.root");

        """.__add__(
            f"""
        size_t max_n_jets = {max_n_jets};
        std::map<int, std::vector<ROOT::RVecI>> permutations = get_permutations_dict(max_n_jets);

        #endif
        """
        )
    )

def define_features(df: ROOT.RDataFrame) -> ROOT.RDataFrame:
    return df.Define(
        "features",
        """
        eval_features(
            permutations.at( std::min(Jet_pt_cut.size(), max_n_jets) ),
            Jet_pt_cut,
            Jet_eta_cut,
            Jet_phi_cut,
            Jet_mass_cut,
            Jet_btagCSVV2_cut,
            Jet_qgl[Jet_mask],
            Electron_pt[Electron_mask],
            Electron_eta[Electron_mask],
            Electron_phi[Electron_mask],
            Electron_mass[Electron_mask],
            Muon_pt[Muon_mask],
            Muon_eta[Muon_mask],
            Muon_phi[Muon_mask],
            Muon_mass[Muon_mask]
        )
        """,
    )

def predict_proba(df: ROOT.RDataFrame) -> ROOT.RDataFrame:
    """get probability scores for every permutation in event"""

    # in inference, odd model applied to even events, while even model to odd events
    # read more about inference in dedicated part of AGC documentation:
    # https://agc.readthedocs.io/en/latest/taskbackground.html#machine-learning-component

    return df.Define(
        "proba",
        """
        bool is_even = (event % 2 == 0);
        const auto& model = (is_even) ? fodd : feven;
        return inference(features, model);
        """,
    )

def infer_output_ml_features(df: ROOT.RDataFrame) -> ROOT.RDataFrame:
    """
    Choose for each feature the best candidate with the highest probability score.
    Results are features of the best four-jet permutation at each event
    """

    df = predict_proba(df)
    for i in range(len(ml_features_config)):
        df = df.Define(f"results{i}", f"features[{i}][ArgMax(proba)]")
    return df

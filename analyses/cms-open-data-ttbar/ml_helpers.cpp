#ifndef ML_HELPERS
#define ML_HELPERS

#include "helpers.h"

#include <cmath>
#include <map>
#include <algorithm>
#include <string>
#include <vector>

#include <TError.h>
#include <ROOT/RVec.hxx>
#include <Math/Vector4D.h>
#include <TMVA/RBDT.hxx>

// copying jet_labels because we need to modify it
std::map<std::string, std::vector<int>> get_permutations (std::string jet_labels) {

    /*
    Builds unique jets permutations.

    jet_labels - an  array of characters, each of them representing jet's label.

    Two w-labels (w-jets), one h-label (hadronic b-tagged jet), one l-label (leptonic b-tagged jet).

    If there is more than 4 jets, other jets gets o-label (stends for other).
    Read more in the official ML part of AGC documentation - https://agc.readthedocs.io/en/latest/taskbackground.html#machine-learning-component

    There is no limitation on characters that are used to denote leptonic or hadronic jets, but W-jets and other jets must be denoted as "w" and "o".
    Leptonic and hadronic jets must be denoted with unique lables. I used "l" and "h" letters as labels.

    Returns indexes of jets of every permutation.
    */

    std::sort(jet_labels.begin(), jet_labels.end());

    /*
    All permutations of indexes are stored in variable permutations.

    Keys of the map are presented as four unique labels: w1 (first W-jet),
                                                         w2 (second W-jet),
                                                         l (leptonic jet),
                                                         h (hadronic jet).

    Every element of the array of indexes corresponds to a unique permutation.
    So the vector {permutations["w1"][i],
                   permutations["w2"][i],
                   permutations["l"][i],
                   permutations["h"][i]} is a i-th permutation of indexes.

    Note:   this function does not produce permutations in the same order as in the reference implementation
    which results in slightly mismatching final histograms.
    possibly because sometimes the order of the 2 w-jets is inverted w.r.t. coffea (ref. implementation)
    */
    std::map<std::string, std::vector<int>> permutations;
    int count = 0, N = jet_labels.size();
    do {
        for (int idx = 0; idx < N; ++idx) { // iterates over all labels of given permutation
            std::string label = std::string(1, jet_labels[idx]); // get labels' character as string
            if (label == "o") continue; // other jets are not stored
            if (label == "w") label+=std::to_string(++count); // gets "w1" or "w2" labels
            permutations[label].push_back(idx); // stores indexes of given permutation
        }
        count = 0; // needs to be reset after every iteration over labels
    } while (std::next_permutation(jet_labels.begin(), jet_labels.end()));


    return permutations;
}

std::map<int, std::vector<ROOT::RVecI>> get_permutations_dict (size_t max_n_jets) {

    // get permutations for different number of jets
    // every event have different number of jets: from 4 to max_n_jets

    std::map<int, std::vector<ROOT::RVecI> > permutations_dict; // number of jets to permutations of jets indexes.
    std::string base = "wwhl";
    for (std::size_t N = 4; N <= max_n_jets; ++N) {
        std::string jet_labels = base + std::string(N-4, 'o');
        std::map<std::string, std::vector<int>> permutations = get_permutations (jet_labels);
        permutations_dict[N] = std::vector<ROOT::RVecI>{permutations["w1"], permutations["w2"], permutations["h"], permutations["l"]};
    }
    return permutations_dict;
}

ROOT::RVec<ROOT::RVecD> eval_features (
    const ROOT::RVec<ROOT::RVecI>& permut_indexes,
    const ROOT::RVecD& jet_pt,
    const ROOT::RVecD& jet_eta,
    const ROOT::RVecD& jet_phi,
    const ROOT::RVecD& jet_mass,
    const ROOT::RVecD& jet_btag,
    const ROOT::RVecD& jet_qgl,
    const ROOT::RVecD& el_pt,
    const ROOT::RVecD& el_eta,
    const ROOT::RVecD& el_phi,
    const ROOT::RVecD& el_mass,
    const ROOT::RVecD& mu_pt,
    const ROOT::RVecD& mu_eta,
    const ROOT::RVecD& mu_phi,
    const ROOT::RVecD& mu_mass
)
{
    using namespace ROOT::VecOps;

    // Part 1. Evaluate Leptons eta, phi, four-momenta
    auto lep_pt = Concatenate(el_pt, mu_pt);
    auto lep_eta = Concatenate(el_eta, mu_eta);
    auto lep_phi = Concatenate(el_phi, mu_phi);
    auto lep_mass = Concatenate(el_mass, mu_mass);
    auto lep_p4 = ConstructP4(lep_pt, lep_eta, lep_phi, lep_mass);

    // Part 2. Evaluating features

    // get jets indexes of each label
    // w1, w2, bH, bL are jets labels

    auto w1_idx = permut_indexes.at(0);
    auto w2_idx = permut_indexes.at(1);
    auto bH_idx = permut_indexes.at(2);
    auto bL_idx = permut_indexes.at(3);


    // Apply indexes to jets.

    auto w1_phi = Take(jet_phi, w1_idx);
    auto w2_phi = Take(jet_phi, w2_idx);
    auto bH_phi = Take(jet_phi, bH_idx);
    auto bL_phi = Take(jet_phi, bL_idx);

    auto w1_eta = Take(jet_eta, w1_idx);
    auto w2_eta = Take(jet_eta, w2_idx);
    auto bH_eta = Take(jet_eta, bH_idx);
    auto bL_eta = Take(jet_eta, bL_idx);

    auto w1_pt = Take(jet_pt, w1_idx); // f8
    auto w2_pt = Take(jet_pt, w2_idx); // f9
    auto bH_pt = Take(jet_pt, bH_idx); // f10
    auto bL_pt = Take(jet_pt, bL_idx); // f11

    auto w1_btag = Take(jet_btag, w1_idx); // f12
    auto w2_btag = Take(jet_btag, w2_idx); // f13
    auto bH_btag = Take(jet_btag, bH_idx); // f14
    auto bL_btag = Take(jet_btag, bL_idx); // f15


    auto w1_qgl = Take(jet_qgl, w1_idx); // f16
    auto w2_qgl = Take(jet_qgl, w2_idx); // f17
    auto bH_qgl = Take(jet_qgl, bH_idx); // f18
    auto bL_qgl = Take(jet_qgl, bL_idx); // f19

    //  jets' four-momenta
    auto w1_p4 = ConstructP4(w1_pt, w1_eta, w1_phi, Take(jet_mass, w1_idx));
    auto w2_p4 = ConstructP4(w2_pt, w2_eta, w2_phi, Take(jet_mass, w2_idx));
    auto bH_p4 = ConstructP4(bH_pt, bH_eta, bH_phi, Take(jet_mass, bH_idx));
    auto bL_p4 = ConstructP4(bL_pt, bL_eta, bL_phi, Take(jet_mass, bL_idx));

    // delta R

    auto dR_lepton_bL = sqrt(pow(lep_eta.at(0)-bL_eta,2.)+pow(lep_phi.at(0)-bL_phi,2.)); // f0
    auto dR_w1w2 = sqrt(pow(w1_eta-w2_eta,2.)+pow(w1_phi-w2_phi,2.)); // f1
    auto dR_w1bH = sqrt(pow(w1_eta-bH_eta,2.)+pow(w1_phi-bH_phi,2.)); // f2
    auto dR_w2bH = sqrt(pow(w2_eta-bH_eta,2.)+pow(w2_phi-bH_phi,2.)); // f3

    // combined mass
    auto M_lepton_bL = Map(bL_p4+lep_p4.at(0), [] (const ROOT::Math::PxPyPzMVector &p) {return p.M();}); //f4
    auto M_w1w2 = Map(w1_p4+w2_p4, [] (const ROOT::Math::PxPyPzMVector &p) {return p.M();}); //f5
    auto M_w1w2bH = Map(w1_p4+w2_p4+bH_p4, [] (const ROOT::Math::PxPyPzMVector &p) {return p.M();}); //f6

    // combined transverse momentum
    auto Pt_w1w2bH = Map(w1_p4+w2_p4+bH_p4, [] (const ROOT::Math::PxPyPzMVector &p) {return p.Pt();}); //f7

    // put all features into single vector
    return ROOT::RVec<ROOT::RVecD> ({
        dR_lepton_bL, dR_w1w2, dR_w1bH, dR_w2bH, M_lepton_bL, M_w1w2, M_w1w2bH, Pt_w1w2bH,
        w1_pt, w2_pt, bH_pt, bL_pt, w1_btag, w2_btag, bH_btag, bL_btag, w1_qgl, w2_qgl, bH_qgl, bL_qgl
    });

}

ROOT::RVecF inference(const ROOT::RVec<ROOT::RVecD> &features, const TMVA::Experimental::RBDT &bdt) {

    size_t npermutations = features.at(0).size();
    ROOT::RVecF res(npermutations);
    size_t nfeatures = features.size();
    ROOT::RVecF input(nfeatures);

    for (std::size_t i = 0; i < npermutations; ++i) {
        for (std::size_t j = 0; j < nfeatures; ++j) {
            input[j] = features.at(j).at(i);
        }
        res[i] = bdt.Compute(input)[0];
    }

    return res;
}

#endif // ML_HELPERS

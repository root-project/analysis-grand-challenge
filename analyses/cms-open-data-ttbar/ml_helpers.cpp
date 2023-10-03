#include "fastforest.h"
#include <cmath>
#include <assert.h>
#include <map>
#include <algorithm>
#include "ROOT/RVec.hxx"

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
        count = 0; // needs to be reset after every itaration over labels
    } while (std::next_permutation(jet_labels.begin(), jet_labels.end()));
    return permutations;
}

std::map<int, std::vector<ROOT::RVecI>> get_permutations_dict (size_t max_n_jets) {

    // get permutations for different number of jets
    // every event have different number of jets: from 4 to max_n_jets
    
    std::map<int, std::vector<ROOT::RVecI> > permutations_dict; // number of jets to permutations of jets indexes.
    std::string base = "wwhl";
    for (int N = 4; N <= max_n_jets; ++N) {
        std::string jet_labels = base + std::string(N-4, 'o'); 
        std::map<std::string, std::vector<int>> permutations = get_permutations (jet_labels);
        permutations_dict[N] = std::vector<ROOT::RVecI>{permutations["w1"], permutations["w2"], permutations["h"], permutations["l"]};
    }
    return permutations_dict;
}

std::map<std::string, fastforest::FastForest> get_fastforests (const std::string& path_to_models, size_t nfeatures) {

    std::vector<std::string> feature_names(nfeatures);
    for (int i = 0; i < nfeatures; ++i) {
        feature_names[i] = "f"+std::to_string(i);
    }

    auto fodd = fastforest::load_txt(path_to_models+"odd.txt", feature_names);
    auto feven = fastforest::load_txt(path_to_models+"even.txt", feature_names);
    return {{"even",feven}, {"odd", fodd}};
}




    
ROOT::RVecF inference(const ROOT::VecOps::RVec<ROOT::RVecD> &features, const fastforest::FastForest &forest, bool check_features=false) {

    size_t npermutations = features.at(0).size();
    size_t nfeatures = features.size();
    ROOT::RVecF res(npermutations);
    float input[nfeatures];

    if (check_features) {
        for (int i = 0; i < nfeatures; ++i) {
            assert(features.at(i).size() == npermutations);
        }
    }
    
    for (int i = 0; i < npermutations; ++i) {
        for (int j = 0; j < nfeatures; ++j) {
            input[j] = features.at(j).at(i);
        }
        float score = forest(input, 0.0F);
        res[i] = 1./(1.+std::exp(-score));
    }

    return res;
}

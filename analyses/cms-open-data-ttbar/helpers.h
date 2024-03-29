#ifndef HELPERS
#define HELPERS

#include <random>
#include <string>
#include <algorithm>

#include "ROOT/RVec.hxx"
#include "TH1D.h"
#include "TRandom3.h"
#include <Math/Vector4D.h>

// functions creating systematic variations
inline double random_gaus()
{
   thread_local std::random_device rd{};
   thread_local std::mt19937 gen{rd()};
   thread_local std::normal_distribution<double> d{1, 0.05};
   return d(gen);
}

inline ROOT::RVecF jet_pt_resolution(std::size_t size)
{
   // normal distribution with 5% variations, shape matches jets
   ROOT::RVecF res(size);
   std::generate(std::begin(res), std::end(res), []()
                 { return random_gaus(); });
   return res;
}

inline float pt_scale_up()
{
   return 1.03;
}

inline ROOT::RVecF btag_weight_variation(const ROOT::RVecF &jet_pt)
{
   // weight variation depending on i-th jet pT (7.5% as default value, multiplied by i-th jet pT / 50 GeV)
   ROOT::RVecF res;
   res.reserve(std::min(4ul, jet_pt.size()));
   for (float pt : Take(jet_pt, 4))
   {
      res.push_back(1 + .075 * pt / 50);
      res.push_back(1 - .075 * pt / 50);
   }
   return res;
}

inline ROOT::RVecF flat_variation()
{
   // 2.5% weight variations
   return 1 + ROOT::RVecF({.025, -.025});
}

inline ROOT::RVec<ROOT::Math::PxPyPzMVector> ConstructP4 (const ROOT::RVecD & Pt, const ROOT::RVecD & Eta, const ROOT::RVecD & Phi, const ROOT::RVecD & M)
{

   return ROOT::VecOps::Construct<ROOT::Math::PxPyPzMVector>(
                ROOT::VecOps::Construct<ROOT::Math::PtEtaPhiMVector>(
                    Pt,
                    Eta,
                    Phi,
                    M
                )
          );
}

#endif

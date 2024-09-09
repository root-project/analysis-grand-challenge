#ifndef HELPERS
#define HELPERS

#include <random>
#include <string>
#include <algorithm>

#include "ROOT/RVec.hxx"
#include "TH1D.h"
#include "TRandom3.h"
#include <Math/Vector4D.h>

inline ROOT::RVecF jet_pt_resolution(const ROOT::RVecF &jet_pt)
{
   // normal distribution with 5% variations, shape matches jets
   ROOT::RVecF res(jet_pt.size());
   // We need to create some pseudo-randomness, it should be thread-safe and at the same time do not depend on RNG. We use the fact that [jet_pt is in GeV....]. 
   // We then use the gaussian quantile to compute the resolution according to the input mean and sigma, using the random bits from the floating-point values.
   double mean = 1.;
   double sigma = 0.05;
   for (std::size_t i = 0; i < jet_pt.size(); ++i) {
      res[i] = mean + ROOT::Math::gaussian_quantile(static_cast<double>(0.001 * (static_cast<int>(jet_pt[i] * 1000) % 1000)) + 0.0005, sigma);
   }

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

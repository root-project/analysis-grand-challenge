#include <string>

#include "ROOT/RVec.hxx"
#include "TH1D.h"
#include "TRandom3.h"

// functions creating systematic variations
inline TRandom &get_thread_local_trandom()
{
   thread_local TRandom rng;
   rng.SetSeed(gRandom->Integer(1000));
   return rng;
}

ROOT::RVecF jet_pt_resolution(std::size_t size)
{
   // normal distribution with 5% variations, shape matches jets
   ROOT::RVecF res(size);
   std::generate(std::begin(res), std::end(res), []()
                 { return get_thread_local_trandom().Gaus(1, 0.05); });
   return res;
}

float pt_scale_up()
{
   return 1.03;
}

ROOT::RVecF btag_weight_variation(const ROOT::RVecF &jet_pt)
{
   // weight variation depending on i-th jet pT (7.5% as default value, multiplied by i-th jet pT / 50 GeV)
   ROOT::RVecF res;
   for (const float &pt : ROOT::VecOps::Take(jet_pt, 4))
   {
      res.push_back(1 + .075 * pt / 50);
      res.push_back(1 - .075 * pt / 50);
   }
   return res;
}

ROOT::RVecF flat_variation()
{
   // 2.5% weight variations
   return 1 + ROOT::RVecF({.025, -.025});
}

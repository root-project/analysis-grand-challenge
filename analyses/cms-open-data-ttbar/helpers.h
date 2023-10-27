#include <random>
#include <string>

#include "ROOT/RVec.hxx"
#include "TH1D.h"
#include "TRandom3.h"
#include <Math/Vector4D.h>


#ifndef HELPERS

// functions creating systematic variations
double random_gaus()
{
   thread_local std::random_device rd{};
   thread_local std::mt19937 gen{rd()};
   thread_local std::normal_distribution<double> d{1, 0.05};
   return d(gen);
}

ROOT::RVecF jet_pt_resolution(std::size_t size)
{
   // normal distribution with 5% variations, shape matches jets
   ROOT::RVecF res(size);
   std::generate(std::begin(res), std::end(res), []()
                 { return random_gaus(); });
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

ROOT::RVec<ROOT::Math::PxPyPzMVector> ConstructP4 (const ROOT::RVecD & Pt, const ROOT::RVecD & Eta, const ROOT::RVecD & Phi, const ROOT::RVecD & M)
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

// to protect from redeclaration error
#define HELPERS

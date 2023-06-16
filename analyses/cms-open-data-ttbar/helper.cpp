#include <string>

#include "ROOT/RVec.hxx"
#include "TH1D.h"
#include "TRandom3.h"

// functions slicing histograms

// accept numbers of bins as parameters
TH1D SliceHisto(const TH1D &h, int xfirst, int xlast)
{

   // do slice in xfirst:xlast including xfirst and xlast
   TH1D res((std::string("h_sliced_") + h.GetTitle()).c_str(), h.GetTitle(), xlast - xfirst,
            h.GetXaxis()->GetBinLowEdge(xfirst), h.GetXaxis()->GetBinUpEdge(xlast - 1));
   // note that histogram arrays are : [ undeflow, bin1, bin2,....., binN, overflow]
   std::copy(h.GetArray() + xfirst, h.GetArray() + xlast, res.GetArray() + 1);
   // set correct underflow/overflows
   res.SetBinContent(0, h.Integral(0, xfirst - 1));                              // set underflow value
   res.SetBinContent(res.GetNbinsX() + 1, h.Integral(xlast, h.GetNbinsX() + 1)); // set overflow value

   return res;
}

// accept axis limits as parameters
TH1D Slice(TH1D &h, double low_edge, double high_edge)
{
   int xfirst = h.FindBin(low_edge);
   int xlast = h.FindBin(high_edge);
   return SliceHisto(h, xfirst, xlast);
}

// functions creating systematic variations
inline TRandom &get_thread_local_trandom()
{
   thread_local TRandom rng;
   rng.SetSeed(gRandom->Integer(1000));
   return rng;
}

ROOT::VecOps::RVec<float> jet_pt_resolution(std::size_t size)
{
   // normal distribution with 5% variations, shape matches jets
   ROOT::VecOps::RVec<float> res(size);
   std::generate(std::begin(res), std::end(res), []()
                 { return get_thread_local_trandom().Gaus(1, 0.05); });
   return res;
}

float pt_scale_up()
{
   return 1.03;
}

ROOT::VecOps::RVec<float> btag_weight_variation(const ROOT::VecOps::RVec<float> &jet_pt)
{
   // weight variation depending on i-th jet pT (7.5% as default value, multiplied by i-th jet pT / 50 GeV)
   ROOT::VecOps::RVec<float> res;
   for (const float &pt : ROOT::VecOps::Take(jet_pt, 4))
   {
      res.push_back(1 + .075 * pt / 50);
      res.push_back(1 - .075 * pt / 50);
   }
   return res;
}

ROOT::VecOps::RVec<float> flat_variation()
{
   // 2.5% weight variations
   return 1 + ROOT::VecOps::RVec<float>({.025, -.025});
}

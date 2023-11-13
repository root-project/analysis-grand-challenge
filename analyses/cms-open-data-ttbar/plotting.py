import ROOT
from ml import ml_features_config
from utils import AGCResult


def save_plots(results: list[AGCResult]):
    width = 2160
    height = 2160
    c = ROOT.TCanvas("c", "c", width, height)
    ROOT.gStyle.SetPalette(ROOT.kRainBow)

    # Region 1 stack
    hlist = [r.histo for r in results if r.region == "4j1b" and r.variation == "nominal"]
    hlist = [h.Clone().Rebin(2) for h in hlist]
    hs = ROOT.THStack("j4b1", ">=4 jets, 1 b-tag; H_{T} [GeV]")
    for h in hlist:
        hs.Add(h)
    hs.Draw("hist pfc plc")
    c.Draw()
    x = hs.GetXaxis()
    x.SetRangeUser(120, x.GetXmax())
    x.SetTitleOffset(1.5)
    x.CenterTitle()
    c.BuildLegend(0.65, 0.7, 0.9, 0.9)
    c.SaveAs("reg1.png")

    # Region 2 stack
    hlist = [r.histo for r in results if r.region == "4j2b" and r.variation == "nominal"]
    hs = ROOT.THStack("j4b2", ">=4 jets, 2 b-tag; m_{bjj} [GeV]")
    for h in hlist:
        hs.Add(h)
    hs.Draw("hist pfc plc")
    c.Draw()
    x = hs.GetXaxis()
    x.SetTitleOffset(1.5)
    x.CenterTitle()
    c.BuildLegend(0.65, 0.7, 0.9, 0.9)
    c.SaveAs("reg2.png")

    # b-tag variations
    btag_variations = [
        "nominal",
        "btag_var_0_up",
        "btag_var_1_up",
        "btag_var_2_up",
        "btag_var_3_up",
    ]

    def btag_filter(r):
        return r.region == "4j1b" and r.process == "ttbar" and r.variation in btag_variations

    hlist = [r.histo for r in results if btag_filter(r)]
    hlist = [h.Clone().Rebin(2) for h in hlist]
    hs = ROOT.THStack("j4b1btag", "btag-variations ; H_{T} [GeV]")
    for h, name in zip(hlist, btag_variations):
        h.SetLineWidth(4)
        h.SetTitle(name)
        hs.Add(h)
    hs.Draw("hist nostack plc")
    c.Draw()
    x = hs.GetXaxis()
    x.SetRangeUser(120, x.GetXmax())
    x.SetTitleOffset(1.5)
    x.CenterTitle()
    c.BuildLegend(0.65, 0.7, 0.9, 0.9)
    c.SaveAs("btag.png")

    # Jet energy variations
    jet_variations = ["nominal", "pt_scale_up", "pt_res_up"]

    def jet_filter(r):
        return r.region == "4j2b" and r.process == "ttbar" and r.variation in jet_variations

    hlist = [r.histo for r in results if jet_filter(r)]
    hs = ROOT.THStack("4j2bjet", "Jet energy variations ; m_{bjj} [GeV]")
    for h, name in zip(hlist, jet_variations):
        h.SetFillColor(0)
        h.SetLineWidth(4)
        h.SetTitle(name)
        hs.Add(h)
    hs.Draw("hist nostack plc")
    c.Draw()
    x = hs.GetXaxis()
    x.SetRangeUser(0, 600)
    x.SetTitleOffset(1.5)
    x.CenterTitle()
    c.BuildLegend(0.65, 0.7, 0.9, 0.9)
    c.SaveAs("jet.png")


def save_ml_plots(results: list[AGCResult]):
    width = 2160
    height = 2160
    c = ROOT.TCanvas("c", "c", width, height)

    for i, feature in enumerate(ml_features_config):
        hlist = [r.histo for r in results if r.variation == "nominal" and r.region == feature]
        hs = ROOT.THStack("features", feature.title)
        for h in hlist:
            hs.Add(h)
        hs.Draw("hist pfc plc")
        c.BuildLegend()
        c.Print("features.pdf" + (i == 0) * "(" + (i + 1 == len(ml_features_config)) * ")")

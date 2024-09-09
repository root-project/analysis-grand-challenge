import os

import ROOT
from statistical_inference.agc_sample import AGCSample
from statistical_inference.drawer import DrawModel, Visualization
from statistical_inference.rebinning_tool import RebinningTool


def fit_histograms(filename=""):
    print(f"Starting fitting process using input file: {filename}")
    
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)

    os.makedirs("statistical_data", exist_ok=True)

    rebinning_tool = RebinningTool()  # Rebinning as extra step
    rebinning_tool.set_xmin(110)
    rebinning_tool.set_rebin(2)
    rebinning_tool.set_input_path(filename)
    rebinning_tool.set_output_path("statistical_data/temp_histos.root")
    rebinning_tool.apply_rebinning()

    file = ROOT.TFile("statistical_data/HistFactoryExtra.root", "RECREATE")  # cleaning the file
    file.Close()  # create separate file for each file -> two much files

    input_file = "statistical_data/temp_histos.root"  # input file for futher defined histograms

    meas = ROOT.RooStats.HistFactory.Measurement("meas", "meas")
    meas.SetLumi(1.0)
    meas.SetLumiRelErr(0.0)

    # ! main difference
    lumi_systematics = (ROOT.RooStats.HistFactory.OverallSys())  # we define it one time -> dont need to redefine it
    lumi_systematics.SetName("Lumi")
    lumi_systematics.SetLow(0.97)
    lumi_systematics.SetHigh(1.03)

    channel = ROOT.RooStats.HistFactory.Channel("channel_4j1b_CR")
    channel.SetData("4j1b_pseudodata", input_file)
    channel.SetStatErrorConfig(0.001, "Gaussian")

    ttbar = AGCSample("ttbar", "4j1b_ttbar", input_file)
    ttbar.SetSystematicsInputFile(input_file)  # if use SetInputFile -> change input file for nominal histogram
    # we can put input file from sample data definition as default

    ttbar.AddOverallSys(lumi_systematics)
    ttbar.ActivateStatError()  # ROOT does not have DisableStatError() method, so we have to keep this :(
    # "histogram_down" can be skipped here:
    ttbar.AddNormPlusShapeHistoSys("ME_variation", histoname_up="4j1b_ttbar_ME_var")
    ttbar.AddNormPlusShapeHistoSys(
        "PS_variation", histoname_up="4j1b_ttbar_PS_var"
    )  # histogram path is specified as "" in default case
    ttbar.AddNormPlusShapeHistoSys(
        "tt_scale_var",
        histoname_up="4j1b_ttbar_scaleup",
        histoname_down="4j1b_ttbar_scaledown",
    )
    ttbar.AddNormPlusShapeHistoSys(
        "jet_energy_scale", histoname_up="4j1b_ttbar_pt_scale_up"
    )
    ttbar.AddNormPlusShapeHistoSys(
        "jet_energy_res", histoname_up="4j1b_ttbar_pt_res_up"
    )
    ttbar.AddNormPlusShapeHistoSys(
        "b_tag_NP_1",
        histoname_up="4j1b_ttbar_btag_var_0_up",
        histoname_down="4j1b_ttbar_btag_var_0_down",
    )
    ttbar.AddNormPlusShapeHistoSys(
        "b_tag_NP_2",
        histoname_up="4j1b_ttbar_btag_var_1_up",
        histoname_down="4j1b_ttbar_btag_var_1_down",
    )
    ttbar.AddNormPlusShapeHistoSys(
        "b_tag_NP_3",
        histoname_up="4j1b_ttbar_btag_var_2_up",
        histoname_down="4j1b_ttbar_btag_var_2_down",
    )
    ttbar.AddNormPlusShapeHistoSys(
        "b_tag_NP_4",
        histoname_up="4j1b_ttbar_btag_var_3_up",
        histoname_down="4j1b_ttbar_btag_var_3_down",
    )
    ttbar.AddNormFactor("ttbar_norm", 1, 0, 10)
    channel.AddSample(ttbar.GetHistFactorySample())

    wjets = AGCSample("wjets", "4j1b_wjets", input_file)
    wjets.SetSystematicsInputFile(input_file)
    wjets.ActivateStatError()
    wjets.AddOverallSys(lumi_systematics)
    wjets.AddNormPlusShapeHistoSys(
        "jet_energy_scale", histoname_up="4j1b_wjets_pt_scale_up"
    )
    wjets.AddNormPlusShapeHistoSys(
        "jet_energy_res", histoname_up="4j1b_wjets_pt_res_up"
    )
    wjets.AddNormPlusShapeHistoSys(
        "b_tag_NP_1",
        histoname_up="4j1b_wjets_btag_var_0_up",
        histoname_down="4j1b_wjets_btag_var_0_down",
    )
    wjets.AddNormPlusShapeHistoSys(
        "b_tag_NP_2",
        histoname_up="4j1b_wjets_btag_var_1_up",
        histoname_down="4j1b_wjets_btag_var_1_down",
    )
    wjets.AddNormPlusShapeHistoSys(
        "b_tag_NP_3",
        histoname_up="4j1b_wjets_btag_var_2_up",
        histoname_down="4j1b_wjets_btag_var_2_down",
    )
    wjets.AddNormPlusShapeHistoSys(
        "b_tag_NP_4",
        histoname_up="4j1b_wjets_btag_var_3_up",
        histoname_down="4j1b_wjets_btag_var_3_down",
    )
    wjets.AddNormPlusShapeHistoSys(
        "w_plus_jets_scale_var",
        histoname_up="4j1b_wjets_scale_var_up",
        histoname_down="4j1b_wjets_scale_var_down",
    )
    channel.AddSample(wjets.GetHistFactorySample())

    single_top_s_chan = AGCSample(
        "single_top_s", "4j1b_single_top_s_chan", input_file
    )
    single_top_s_chan.SetSystematicsInputFile(input_file)
    single_top_s_chan.ActivateStatError()
    single_top_s_chan.AddOverallSys(lumi_systematics)
    single_top_s_chan.AddNormPlusShapeHistoSys(
        "jet_energy_scale", histoname_up="4j1b_single_top_s_chan_pt_scale_up"
    )
    single_top_s_chan.AddNormPlusShapeHistoSys(
        "jet_energy_res", histoname_up="4j1b_single_top_s_chan_pt_res_up"
    )
    single_top_s_chan.AddNormPlusShapeHistoSys(
        "b_tag_NP_1",
        histoname_up="4j1b_single_top_s_chan_btag_var_0_up",
        histoname_down="4j1b_single_top_s_chan_btag_var_0_down",
    )
    single_top_s_chan.AddNormPlusShapeHistoSys(
        "b_tag_NP_2",
        histoname_up="4j1b_single_top_s_chan_btag_var_1_up",
        histoname_down="4j1b_single_top_s_chan_btag_var_1_down",
    )
    single_top_s_chan.AddNormPlusShapeHistoSys(
        "b_tag_NP_3",
        histoname_up="4j1b_single_top_s_chan_btag_var_2_up",
        histoname_down="4j1b_single_top_s_chan_btag_var_2_down",
    )
    single_top_s_chan.AddNormPlusShapeHistoSys(
        "b_tag_NP_4",
        histoname_up="4j1b_single_top_s_chan_btag_var_3_up",
        histoname_down="4j1b_single_top_s_chan_btag_var_3_down",
    )

    channel.AddSample(single_top_s_chan.GetHistFactorySample())

    single_top_t_chan = AGCSample(
        "single_top_t", "4j1b_single_top_t_chan", input_file
    )
    single_top_t_chan.SetSystematicsInputFile(input_file)
    single_top_t_chan.ActivateStatError()
    single_top_t_chan.AddOverallSys(lumi_systematics)
    single_top_t_chan.AddNormPlusShapeHistoSys(
        "jet_energy_scale", histoname_up="4j1b_single_top_t_chan_pt_scale_up"
    )
    single_top_t_chan.AddNormPlusShapeHistoSys(
        "jet_energy_res", histoname_up="4j1b_single_top_t_chan_pt_res_up"
    )
    single_top_t_chan.AddNormPlusShapeHistoSys(
        "b_tag_NP_1",
        histoname_up="4j1b_single_top_t_chan_btag_var_0_up",
        histoname_down="4j1b_single_top_t_chan_btag_var_0_down",
    )
    single_top_t_chan.AddNormPlusShapeHistoSys(
        "b_tag_NP_2",
        histoname_up="4j1b_single_top_t_chan_btag_var_1_up",
        histoname_down="4j1b_single_top_t_chan_btag_var_1_down",
    )
    single_top_t_chan.AddNormPlusShapeHistoSys(
        "b_tag_NP_3",
        histoname_up="4j1b_single_top_t_chan_btag_var_2_up",
        histoname_down="4j1b_single_top_t_chan_btag_var_2_down",
    )
    single_top_t_chan.AddNormPlusShapeHistoSys(
        "b_tag_NP_4",
        histoname_up="4j1b_single_top_t_chan_btag_var_3_up",
        histoname_down="4j1b_single_top_t_chan_btag_var_3_down",
    )

    channel.AddSample(single_top_t_chan.GetHistFactorySample())

    single_top_tW = AGCSample(
        "single_top_tW", "4j1b_single_top_tW", input_file
    )
    single_top_tW.SetSystematicsInputFile(input_file)
    single_top_tW.ActivateStatError()
    single_top_tW.AddOverallSys(lumi_systematics)
    single_top_tW.AddNormPlusShapeHistoSys(
        "jet_energy_scale", histoname_up="4j1b_single_top_tW_pt_scale_up"
    )
    single_top_tW.AddNormPlusShapeHistoSys(
        "jet_energy_res", histoname_up="4j1b_single_top_tW_pt_res_up"
    )
    single_top_tW.AddNormPlusShapeHistoSys(
        "b_tag_NP_1",
        histoname_up="4j1b_single_top_tW_btag_var_0_up",
        histoname_down="4j1b_single_top_tW_btag_var_0_down",
    )
    single_top_tW.AddNormPlusShapeHistoSys(
        "b_tag_NP_2",
        histoname_up="4j1b_single_top_tW_btag_var_1_up",
        histoname_down="4j1b_single_top_tW_btag_var_1_down",
    )
    single_top_tW.AddNormPlusShapeHistoSys(
        "b_tag_NP_3",
        histoname_up="4j1b_single_top_tW_btag_var_2_up",
        histoname_down="4j1b_single_top_tW_btag_var_2_down",
    )
    single_top_tW.AddNormPlusShapeHistoSys(
        "b_tag_NP_4",
        histoname_up="4j1b_single_top_tW_btag_var_3_up",
        histoname_down="4j1b_single_top_tW_btag_var_3_down",
    )

    channel.AddSample(single_top_tW.GetHistFactorySample())

    meas.AddChannel(channel)

    channel_2b = ROOT.RooStats.HistFactory.Channel("channel_4j2b_SR")
    channel_2b.SetData("4j2b_pseudodata", input_file)
    channel_2b.SetStatErrorConfig(0.001, "Gaussian")

    ttbar = AGCSample("ttbar", "4j2b_ttbar", input_file)
    ttbar.AddOverallSys(lumi_systematics)
    ttbar.SetSystematicsInputFile(input_file)
    ttbar.ActivateStatError()
    ttbar.AddNormPlusShapeHistoSys(
        "ME_variation", histoname_up="4j2b_ttbar_ME_var"
    )
    ttbar.AddNormPlusShapeHistoSys(
        "PS_variation", histoname_up="4j2b_ttbar_PS_var"
    )
    ttbar.AddNormPlusShapeHistoSys(
        "tt_scale_var",
        histoname_up="4j2b_ttbar_scaleup",
        histoname_down="4j2b_ttbar_scaledown",
    )
    ttbar.AddNormPlusShapeHistoSys(
        "jet_energy_scale", histoname_up="4j2b_ttbar_pt_scale_up"
    )
    ttbar.AddNormPlusShapeHistoSys(
        "jet_energy_res", histoname_up="4j2b_ttbar_pt_res_up"
    )
    ttbar.AddNormPlusShapeHistoSys(
        "b_tag_NP_1",
        histoname_up="4j2b_ttbar_btag_var_0_up",
        histoname_down="4j2b_ttbar_btag_var_0_down",
    )
    ttbar.AddNormPlusShapeHistoSys(
        "b_tag_NP_2",
        histoname_up="4j2b_ttbar_btag_var_1_up",
        histoname_down="4j2b_ttbar_btag_var_1_down",
    )
    ttbar.AddNormPlusShapeHistoSys(
        "b_tag_NP_3",
        histoname_up="4j2b_ttbar_btag_var_2_up",
        histoname_down="4j2b_ttbar_btag_var_2_down",
    )
    ttbar.AddNormPlusShapeHistoSys(
        "b_tag_NP_4",
        histoname_up="4j2b_ttbar_btag_var_3_up",
        histoname_down="4j2b_ttbar_btag_var_3_down",
    )
    ttbar.AddNormFactor("ttbar_norm", 1, 0, 10)
    channel_2b.AddSample(ttbar.GetHistFactorySample())

    wjets = AGCSample("wjets", "4j2b_wjets", input_file)
    wjets.SetSystematicsInputFile(input_file)
    wjets.ActivateStatError()
    wjets.AddOverallSys(lumi_systematics)
    wjets.AddNormPlusShapeHistoSys(
        "jet_energy_scale", histoname_up="4j2b_wjets_pt_scale_up"
    )
    wjets.AddNormPlusShapeHistoSys(
        "jet_energy_res", histoname_up="4j2b_wjets_pt_res_up"
    )
    wjets.AddNormPlusShapeHistoSys(
        "b_tag_NP_1",
        histoname_up="4j2b_wjets_btag_var_0_up",
        histoname_down="4j2b_wjets_btag_var_0_down",
    )
    wjets.AddNormPlusShapeHistoSys(
        "b_tag_NP_2",
        histoname_up="4j2b_wjets_btag_var_1_up",
        histoname_down="4j2b_wjets_btag_var_1_down",
    )
    wjets.AddNormPlusShapeHistoSys(
        "b_tag_NP_3",
        histoname_up="4j2b_wjets_btag_var_2_up",
        histoname_down="4j2b_wjets_btag_var_2_down",
    )
    wjets.AddNormPlusShapeHistoSys(
        "b_tag_NP_4",
        histoname_up="4j2b_wjets_btag_var_3_up",
        histoname_down="4j2b_wjets_btag_var_3_down",
    )
    wjets.AddNormPlusShapeHistoSys(
        "w_plus_jets_scale_var",
        histoname_up="4j2b_wjets_scale_var_up",
        histoname_down="4j2b_wjets_scale_var_down",
    )
    channel_2b.AddSample(wjets.GetHistFactorySample())

    single_top_s_chan = AGCSample(
        "single_top_s", "4j2b_single_top_s_chan", input_file
    )
    single_top_s_chan.SetSystematicsInputFile(input_file)
    single_top_s_chan.ActivateStatError()
    single_top_s_chan.AddOverallSys(lumi_systematics)
    single_top_s_chan.AddNormPlusShapeHistoSys(
        "jet_energy_scale", histoname_up="4j2b_single_top_s_chan_pt_scale_up"
    )
    single_top_s_chan.AddNormPlusShapeHistoSys(
        "jet_energy_res", histoname_up="4j2b_single_top_s_chan_pt_res_up"
    )
    single_top_s_chan.AddNormPlusShapeHistoSys(
        "b_tag_NP_1",
        histoname_up="4j2b_single_top_s_chan_btag_var_0_up",
        histoname_down="4j2b_single_top_s_chan_btag_var_0_down",
    )
    single_top_s_chan.AddNormPlusShapeHistoSys(
        "b_tag_NP_2",
        histoname_up="4j2b_single_top_s_chan_btag_var_1_up",
        histoname_down="4j2b_single_top_s_chan_btag_var_1_down",
    )
    single_top_s_chan.AddNormPlusShapeHistoSys(
        "b_tag_NP_3",
        histoname_up="4j2b_single_top_s_chan_btag_var_2_up",
        histoname_down="4j2b_single_top_s_chan_btag_var_2_down",
    )
    single_top_s_chan.AddNormPlusShapeHistoSys(
        "b_tag_NP_4",
        histoname_up="4j2b_single_top_s_chan_btag_var_3_up",
        histoname_down="4j2b_single_top_s_chan_btag_var_3_down",
    )
    channel_2b.AddSample(single_top_s_chan.GetHistFactorySample())

    single_top_t_chan = AGCSample(
        "single_top_t", "4j2b_single_top_t_chan", input_file
    )
    single_top_t_chan.SetSystematicsInputFile(input_file)
    single_top_t_chan.ActivateStatError()
    single_top_t_chan.AddOverallSys(lumi_systematics)
    single_top_t_chan.AddNormPlusShapeHistoSys(
        "jet_energy_scale", histoname_up="4j2b_single_top_t_chan_pt_scale_up"
    )
    single_top_t_chan.AddNormPlusShapeHistoSys(
        "jet_energy_res", histoname_up="4j2b_single_top_t_chan_pt_res_up"
    )
    single_top_t_chan.AddNormPlusShapeHistoSys(
        "b_tag_NP_1",
        histoname_up="4j2b_single_top_t_chan_btag_var_0_up",
        histoname_down="4j2b_single_top_t_chan_btag_var_0_down",
    )
    single_top_t_chan.AddNormPlusShapeHistoSys(
        "b_tag_NP_2",
        histoname_up="4j2b_single_top_t_chan_btag_var_1_up",
        histoname_down="4j2b_single_top_t_chan_btag_var_1_down",
    )
    single_top_t_chan.AddNormPlusShapeHistoSys(
        "b_tag_NP_3",
        histoname_up="4j2b_single_top_t_chan_btag_var_2_up",
        histoname_down="4j2b_single_top_t_chan_btag_var_2_down",
    )
    single_top_t_chan.AddNormPlusShapeHistoSys(
        "b_tag_NP_4",
        histoname_up="4j2b_single_top_t_chan_btag_var_3_up",
        histoname_down="4j2b_single_top_t_chan_btag_var_3_down",
    )

    channel_2b.AddSample(single_top_t_chan.GetHistFactorySample())

    single_top_tW = AGCSample(
        "single_top_tW", "4j2b_single_top_tW", input_file
    )
    single_top_tW.ActivateStatError()
    single_top_tW.SetSystematicsInputFile(input_file)
    single_top_tW.AddOverallSys(lumi_systematics)
    single_top_tW.AddNormPlusShapeHistoSys(
        "jet_energy_scale", histoname_up="4j2b_single_top_tW_pt_scale_up"
    )
    single_top_tW.AddNormPlusShapeHistoSys(
        "jet_energy_res", histoname_up="4j2b_single_top_tW_pt_res_up"
    )
    single_top_tW.AddNormPlusShapeHistoSys(
        "b_tag_NP_1",
        histoname_up="4j2b_single_top_tW_btag_var_0_up",
        histoname_down="4j2b_single_top_tW_btag_var_0_down",
    )
    single_top_tW.AddNormPlusShapeHistoSys(
        "b_tag_NP_2",
        histoname_up="4j2b_single_top_tW_btag_var_1_up",
        histoname_down="4j2b_single_top_tW_btag_var_1_down",
    )
    single_top_tW.AddNormPlusShapeHistoSys(
        "b_tag_NP_3",
        histoname_up="4j2b_single_top_tW_btag_var_2_up",
        histoname_down="4j2b_single_top_tW_btag_var_2_down",
    )
    single_top_tW.AddNormPlusShapeHistoSys(
        "b_tag_NP_4",
        histoname_up="4j2b_single_top_tW_btag_var_3_up",
        histoname_down="4j2b_single_top_tW_btag_var_3_down",
    )

    channel_2b.AddSample(single_top_tW.GetHistFactorySample())
    meas.AddChannel(channel_2b)

    meas.SetPOI("ttbar_norm")
    meas.CollectHistograms()

    ws = ROOT.RooStats.HistFactory.MakeModelAndMeasurementFast(meas)

    # Retrieve the ModelConfig
    modelConfig = ws.obj("ModelConfig")

    # Extract the PDF and global observables
    pdf = modelConfig.GetPdf()
    globalObservables = ROOT.RooArgSet(modelConfig.GetGlobalObservables())

    # Perform the fit
    result = pdf.fitTo(
        ws["obsData"],
        Save=True,
        PrintLevel=-1,
        GlobalObservables=globalObservables,
        SumW2Error=False
    )

    vis = Visualization()
    vis.CreateAndSaveResultsPicture("fitted_parameters.png", ws)  # Parameters fit result

    vis.DrawCorrelationMatrix("correlation_matrix.png", result)  # Correlation matrix

    md = DrawModel(meas, ws)
    md.Draw(result)

    # save fit results to ROOT file

    with ROOT.TFile.Open("fitResults.root", "RECREATE") as file:
            result.Write("fitResult")


    result.Print()

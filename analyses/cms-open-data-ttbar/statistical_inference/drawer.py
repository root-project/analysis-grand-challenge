import numpy as np
import ROOT


class Visualization:
    def CreateAndSaveResultsPicture(self, filename, fWorkspace):
        names = []
        values = []
        errors = []

        vars = fWorkspace.allVars()
        # for (auto var : vars)
        for var in vars:
            name = var.GetName()
            if "alpha" in name:
                if fWorkspace.var(name).getVal() == 0:
                    continue
                names += [name.split("alpha_")[-1]]
                values += [fWorkspace.var(name).getVal()]
                errors += [fWorkspace.var(name).getError()]

        num = len(names)

        ROOT.gStyle.SetPalette(1)
        self.c1 = ROOT.TCanvas("c1", "c1", 1200, 600)
        self.c1.SetLeftMargin(0.2)

        self.frame = ROOT.TH2F("self.frame", "", 6, -3, 3, num, 0, num)

        self.frame.Draw("")

        self.box = ROOT.TBox(-2, 0, 2, num)
        # external yellow box, which represent 2 sigma deviation
        self.box.SetFillColorAlpha(ROOT.kYellow - 9, 0.8)
        self.box.Draw("same A")

        self.box_internal = ROOT.TBox(-1, 0, 1, num)
        # internal green box, which represent 1 sigma deviation
        self.box_internal.SetFillColorAlpha(ROOT.kGreen, 0.5)
        self.box_internal.Draw("same A")

        self.axis = ROOT.TGaxis(
            -3, num, 3, num, -3, 3, 510, "-"
        )  # left Y axis, which used for naming parameters

        self.frame.GetYaxis().SetTickLength(
            0.0
        )  # set default axis ticks to 0 -> they will be overwritten by manually created ones

        self.graph = ROOT.TGraph()
        self.lines = []
        self.right_ticks = []
        self.left_ticks = []

        self.frame.GetYaxis().SetLabelSize(0.04)

        for i in range(num):
            self.graph.SetPoint(i, values[i], i + 0.5)
            self.frame.GetYaxis().SetBinLabel(i + 1, names[i])

            self.lines += [
                ROOT.TLine(
                    values[i] - errors[i],
                    i + 0.5,
                    values[i] + errors[i],
                    i + 0.5,
                )
            ]  # line which represent parameter error
            self.lines[-1].SetLineColor(ROOT.kBlack)
            self.lines[-1].SetLineWidth(2)
            self.lines[-1].Draw()

            self.left_ticks += [
                ROOT.TLine(-3, i + 0.5, -2.95, i + 0.5)
            ]  # left tick, paralel to parameter error line
            self.left_ticks[-1].SetLineColor(ROOT.kBlack)
            self.left_ticks[-1].SetLineWidth(1)
            self.left_ticks[-1].Draw()

            self.right_ticks += [
                ROOT.TLine(2.95, i + 0.5, 3, i + 0.5)
            ]  # right tick, paralel to parameter error line
            self.right_ticks[-1].SetLineColor(ROOT.kBlack)
            self.right_ticks[-1].SetLineWidth(1)
            self.right_ticks[-1].Draw()

        self.tl = ROOT.TLine(0, 0, 0, num)
        self.tl.SetLineStyle(2)
        self.tl.Draw()

        self.graph.SetMarkerStyle(20)
        self.graph.SetMarkerSize(1)
        self.graph.SetMarkerColor(ROOT.kBlack)

        self.graph.Draw("P same")

        self.frame.SetStats(0)
        self.axis.SetLabelSize(0.0)
        self.axis.Draw()

        ROOT.gPad.RedrawAxis()

        self.c1.SaveAs(filename)
        self.c1.Draw()

    def DrawCorrelationMatrix(self, filename, result):
        final_parameters = result.floatParsFinal()
        corr_matrix_before = result.correlationMatrix()

        number_of_inter_params = 0

        for i in range(len(final_parameters)):
            par = final_parameters.at(i)
            if "gamma" in par.GetName():  # skip parameters for bins stat error
                continue

            number_of_inter_params += 1

        name = "CorrelationMatrix for fit results"

        n = corr_matrix_before.GetNcols()

        self.hh = ROOT.TH2D(
            name,
            name,
            number_of_inter_params,
            0,
            number_of_inter_params,
            number_of_inter_params,
            0,
            number_of_inter_params,
        )

        internal_index = 0
        for i in range(n):

            par = final_parameters.at(i)
            if "gamma" in par.GetName():
                continue

            internal__internal_index = 0
            for j in range(n):
                par = final_parameters.at(j)
                if "gamma" in par.GetName():
                    continue
                self.hh.Fill(
                    internal_index + 0.5,
                    number_of_inter_params - internal__internal_index - 0.5,
                    corr_matrix_before[i][j],
                )
                internal__internal_index += 1

            self.hh.GetXaxis().SetBinLabel(
                internal_index + 1,
                final_parameters[i].GetName().split("alpha_")[-1],
            )
            self.hh.GetYaxis().SetBinLabel(
                number_of_inter_params - internal_index,
                final_parameters[i].GetName().split("alpha_")[-1],
            )
            internal_index += 1

        self.hh.SetMinimum(-1)
        self.hh.SetMaximum(+1)

        self.c = ROOT.TCanvas(
            "self.c",
            "Canvas",
            number_of_inter_params * 100,
            number_of_inter_params * 60,
        )
        self.hh.Draw("COLZ")
        self.hh.SetStats(0)

        ROOT.gStyle.SetPalette(87)
        palette = self.hh.GetListOfFunctions().FindObject("palette")
        if palette:
            palette.SetX1NDC(0.1)  # Adjust palette position
            palette.SetX2NDC(0.3)  # Adjust palette position

        self.c.SaveAs(filename)
        self.c.Draw()


class DrawModel:

    predefined_colors = [
        ROOT.TColor.GetColor(
            "#3F90DA"
        ),  # check choosen colors here: https://arxiv.org/abs/2107.02270
        ROOT.TColor.GetColor("#FFA90E"),
        ROOT.TColor.GetColor("#BD1F01"),
        ROOT.TColor.GetColor("#94A4A2"),
        ROOT.TColor.GetColor("#832DB6"),
        ROOT.TColor.GetColor("#A96B59"),
        ROOT.TColor.GetColor("#E76300"),
        ROOT.TColor.GetColor("#B9AC70"),
        ROOT.TColor.GetColor("#717581"),
        ROOT.TColor.GetColor("#92DADD"),
        ROOT.TColor.GetColor("#D0E1F9"),  # Light Azure
        ROOT.TColor.GetColor("#E3D6F1"),  # Light Violet
        ROOT.TColor.GetColor("#000000"),  # Black (for contrast)
        ROOT.TColor.GetColor("#C0C0C0"),  # Gray
        ROOT.TColor.GetColor("#F4A3A0"),  # Soft Red
        ROOT.TColor.GetColor("#9AB9F5"),  # Soft Blue
        ROOT.TColor.GetColor("#B5E5B0"),  # Soft Green
        ROOT.TColor.GetColor("#F0A1B2"),  # Soft Magenta
        ROOT.TColor.GetColor("#B0D6D5"),  # Soft Cyan
        ROOT.TColor.GetColor("#F3F9A6"),  # Soft Yellow
    ]

    def GetChannelPrefitGraph(self, number):
        return self.hs_stacks[number]

    def GetChannelPostfitGraph(self, number):
        return self.second_hs_stacks[number]

    def __init__(self, meas, ws):
        self.meas = meas
        self.ws = ws
        self.cv_array = (
            []
        )  # save all ROOT Graphic objects as self to be able show picture in main notebook
        self.hs_stacks = []
        self.sample_histograms = []
        self.bias_graphs = []
        self.bias_second_graphs = []
        self.normal_lines = []
        self.second_histos = []
        self.second_hs_stacks = []
        self.afterfit_histograms = []
        self.error_graphs = []
        self.error_pre_graphs = []

    def get_yields(self, variable, observables, pdf, result, prefit=False):

        yields = np.zeros(variable.numBins())
        yields_uncert = np.zeros(variable.numBins())
        sample_values = {}
        for i_bin in range(variable.numBins()):
            variable.setBin(i_bin)
            bin_width = variable.getBinWidth(i_bin)

            if prefit:  # for prefit original parameters values
                fill_array = np.zeros(
                    (
                        pdf.funcList().size(),
                        len(pdf.getParameters(observables)),
                    )
                )
                all_params = pdf.getParameters(
                    observables
                )  # get parameters, which required only for interested pdf
                param_values = (
                    result.floatParsInit()
                )  # original parameters values

                for i_sample in range(pdf.funcList().size()):
                    sample_yield = ROOT.RooProduct(
                        "tmp",
                        "tmp",
                        [pdf.funcList()[i_sample], pdf.coefList()[i_sample]],
                    )  # get sample from pdf

                    yields[i_bin] += (
                        bin_width * sample_yield.getVal()
                    )  # get total bin content per bin
                    if i_sample in sample_values:
                        sample_values[i_sample] += [
                            bin_width * sample_yield.getVal()
                        ]  # get value for each sample per bin
                    else:
                        sample_values[i_sample] = [
                            bin_width * sample_yield.getVal()
                        ]

                    for j_parameter, par in enumerate(all_params):

                        name = par.GetName()

                        original_ind = param_values.index(
                            param_values.find(name)
                        )  # index of parameter from pdf in original parameters array

                        postfit_cen_val = (
                            par.getVal()
                        )  # saving postfit value to restore it after all calculations
                        cen_val = param_values[
                            original_ind
                        ].getVal()  # getting original value for parameter
                        par_err = param_values[
                            original_ind
                        ].getError()  # getting original error for parameter

                        par.setVal(cen_val + par_err)
                        par_upper_variation = sample_yield.getVal(
                            observables
                        )  # get upper variation for sample
                        par.setVal(cen_val - par_err)
                        par_bottom_variation = sample_yield.getVal(
                            observables
                        )  # get lower variation for sample

                        par.setVal(postfit_cen_val)  # restore postfit value

                        fill_array[i_sample, j_parameter] = (
                            par_upper_variation - par_bottom_variation
                        ) / 2

                total_uncertanties_per_variation = np.sum(
                    fill_array, axis=0
                )  # sum by 0 axis to get total uncertainty per variation
                total_uncertainty = np.sum(
                    np.power(total_uncertanties_per_variation, 2), axis=0
                )  # for postfit covariance matrix = identity matrix, so code be bit simplified
                yields_uncert[i_bin] = (
                    np.sqrt(total_uncertainty) * bin_width
                )  # original values is normalized -> multiply by bin width

            else:
                all_pdf_params = pdf.getParameters(observables)

                required_params = [
                    i for i in all_pdf_params if i.getError() > 0.001
                ]  # get only non zero parameters

                fill_array = np.zeros(
                    (pdf.funcList().size(), len(required_params))
                )

                number_of_parameters = len(required_params)
                cov_matrix = result.reducedCovarianceMatrix(
                    required_params
                )  # reduced covariance matrix for only non 0 parameters

                numpy_cov_matrix = np.zeros(
                    (number_of_parameters, number_of_parameters)
                )

                for i_index in range(
                    number_of_parameters
                ):  # fill and normalize numpy covariance matrix
                    for j_index in range(i_index, number_of_parameters):
                        numpy_cov_matrix[i_index, j_index] = cov_matrix[
                            i_index, j_index
                        ] / np.sqrt(
                            cov_matrix[i_index, i_index]
                            * cov_matrix[j_index, j_index]
                        )
                        numpy_cov_matrix[j_index, i_index] = numpy_cov_matrix[
                            i_index, j_index
                        ]

                for i_sample in range(pdf.funcList().size()):
                    sample_yield = ROOT.RooProduct(
                        "tmp",
                        "tmp",
                        [pdf.funcList()[i_sample], pdf.coefList()[i_sample]],
                    )  # get sample from pdf
                    yields[i_bin] += (
                        bin_width * sample_yield.getVal()
                    )  # get total bin content per bin
                    if i_sample in sample_values:
                        sample_values[i_sample] += [
                            bin_width * sample_yield.getVal()
                        ]  # get value for each sample per bin
                    else:
                        sample_values[i_sample] = [
                            bin_width * sample_yield.getVal()
                        ]

                    for j_parameter, par in enumerate(required_params):

                        cen_val = (
                            par.getVal()
                        )  # getting postfit parameter value
                        par_err = (
                            par.getError()
                        )  # getting postfit parameter error

                        par.setVal(cen_val + par_err)
                        par_upper_variation = sample_yield.getVal(
                            observables
                        )  # getting upper variation
                        par.setVal(cen_val - par_err)
                        par_bottom_variation = sample_yield.getVal(
                            observables
                        )  # getting lower variation

                        par.setVal(cen_val)  # restore postfit value
                        sample_yield.getVal(observables)

                        fill_array[i_sample, j_parameter] = (
                            par_upper_variation - par_bottom_variation
                        ) / 2

                total_uncertanties_per_variation = np.sum(
                    fill_array, axis=0
                )  # sum by 0 axis to get total uncertainty per variation

                total_uncertainty = np.dot(
                    total_uncertanties_per_variation,
                    np.dot(numpy_cov_matrix, total_uncertanties_per_variation),
                )  # F * (C * F)

                yields_uncert[i_bin] = np.sqrt(total_uncertainty) * bin_width

        return yields, yields_uncert, sample_values

    def get_yields_no_fit(
        self, variable, observables, pdf, result
    ):  # for extra variables and not fit
        yields = np.zeros(variable.numBins())
        yields_uncert = np.zeros(variable.numBins())
        sample_values = {}
        for i_bin in range(variable.numBins()):
            variable.setBin(i_bin)
            bin_width = variable.getBinWidth(i_bin)

            required_params = pdf.getParameters(
                observables
            )  # get all parameters required for uncertainty calculation

            fit_result_params = (
                result.floatParsFinal()
            )  # postfit parameters_values

            fill_array = np.zeros(
                (pdf.funcList().size(), len(required_params))
            )

            number_of_parameters = len(required_params)

            parameters_indices_map = {}
            parameters_for_cov_matrix = []

            internal_index = 0
            for i_index, par in enumerate(
                required_params
            ):  # loop over new pdf parameters and setting their values to postfit parameters values
                par_index = fit_result_params.index(
                    fit_result_params.find(par.GetName())
                )
                if par_index == -1:
                    continue
                parameters_indices_map[i_index] = (
                    internal_index  # internal index for copying covariance matrix
                )
                internal_index += 1
                parameter_to_be_copied = fit_result_params[par_index]
                par.setVal(parameter_to_be_copied.getVal())
                par.setError(parameter_to_be_copied.getError())
                parameters_for_cov_matrix += [
                    parameter_to_be_copied
                ]  # save parameters to list to get reduced Covariance martix

            cov_matrix = result.reducedCovarianceMatrix(
                parameters_for_cov_matrix
            )

            numpy_cov_matrix = np.zeros(
                (number_of_parameters, number_of_parameters)
            )

            for i_index in range(
                number_of_parameters
            ):  # filling covariance matrix using values from postfit covariance matrix.
                for j_index in range(i_index, number_of_parameters):
                    if (
                        i_index in parameters_indices_map
                        and j_index in parameters_indices_map
                    ):  # if both parameters in row and column existing in postfit results -> copy and normalize
                        numpy_cov_matrix[i_index, j_index] = cov_matrix[
                            parameters_indices_map[i_index],
                            parameters_indices_map[j_index],
                        ] / np.sqrt(
                            cov_matrix[
                                parameters_indices_map[i_index],
                                parameters_indices_map[i_index],
                            ]
                            * cov_matrix[
                                parameters_indices_map[j_index],
                                parameters_indices_map[j_index],
                            ]
                        )
                        numpy_cov_matrix[j_index, i_index] = numpy_cov_matrix[
                            i_index, j_index
                        ]
                    else:  # if at least one parameter was not presented in postfit results
                        numpy_cov_matrix[i_index, j_index] = int(
                            i_index == j_index
                        )

            """
            expected view of covariance matrix after copying:
            var1, var2, var3 was presented in postfit results
            var4, var5 are new parameters (mostly are stat errors for bins)

                    var1    var2    var3    var4    var5
            var1    1       non0    non0    0       0
            var2    non0    1       non0    0       0
            var3    non0    non0    1       0       0
            var4    0       0       0       1       0
            var5    0       0       0       0       1


            """

            for i_sample in range(
                pdf.funcList().size()
            ):  # same as default one
                sample_yield = ROOT.RooProduct(
                    "tmp",
                    "tmp",
                    [pdf.funcList()[i_sample], pdf.coefList()[i_sample]],
                )
                yields[i_bin] += bin_width * sample_yield.getVal()
                if i_sample in sample_values:
                    sample_values[i_sample] += [
                        bin_width * sample_yield.getVal()
                    ]
                else:
                    sample_values[i_sample] = [
                        bin_width * sample_yield.getVal()
                    ]

                for j_parameter, par in enumerate(required_params):

                    cen_val = par.getVal()
                    par_err = par.getError()

                    par.setVal(cen_val + par_err)
                    par_upper_variation = sample_yield.getVal(observables)
                    par.setVal(cen_val - par_err)
                    par_bottom_variation = sample_yield.getVal(observables)

                    par.setVal(cen_val)
                    sample_yield.getVal(observables)

                    fill_array[i_sample, j_parameter] = (
                        par_upper_variation - par_bottom_variation
                    ) / 2

            total_uncertanties_per_variation = np.sum(fill_array, axis=0)

            total_uncertainty = np.dot(
                total_uncertanties_per_variation,
                np.dot(numpy_cov_matrix, total_uncertanties_per_variation),
            )

            yields_uncert[i_bin] = np.sqrt(total_uncertainty) * bin_width

        return yields, yields_uncert, sample_values

    def Draw(self, result, no_fit=False):
        self.boxes = []
        self.error_boxes = []
        self.error_boxes_prefit = []
        self.data_plots = []
        for channel in self.meas.GetChannels():
            channel_pdf = self.ws[str(channel.GetName()) + "_model"]
            obs_var = self.ws["obs_x_" + str(channel.GetName())]
            observables = ROOT.RooArgSet(obs_var)

            divide_value = 0.3  # divide canvas into 4 pads: prefit histogram, postfit histogram, prefit data/bin_value, postfit data/bin_value

            self.cv_array += [
                ROOT.TCanvas(
                    "canvas" + str(channel.GetName()),
                    "canvas" + str(channel.GetName()),
                    1500,
                    600,
                )
            ]
            pad1_upper = ROOT.TPad(
                "pad1_upper" + str(channel.GetName()),
                "pad1_upper" + str(channel.GetName()),
                0,
                divide_value,
                0.5,
                1,
            )
            pad1_upper.Draw()
            pad1_bottom = ROOT.TPad(
                "pad1_bottom" + str(channel.GetName()),
                "pad1_bottom" + str(channel.GetName()),
                0,
                0,
                0.5,
                divide_value,
            )
            pad1_bottom.Draw()

            pad2_upper = ROOT.TPad(
                "pad2_upper" + str(channel.GetName()),
                "pad2_upper" + str(channel.GetName()),
                0.5,
                divide_value,
                1,
                1,
            )
            pad2_upper.Draw()
            pad2_bottom = ROOT.TPad(
                "pad2_bottom" + str(channel.GetName()),
                "pad2_bottom" + str(channel.GetName()),
                0.5,
                0,
                1,
                divide_value,
            )
            pad2_bottom.Draw()
            pad1_upper.cd()

            prefit_yields, prefit_unc, prefit_sample_values = self.get_yields(
                obs_var, observables, channel_pdf, result, prefit=True
            )  # get prefit uncert for histogram

            self.hs_stacks += [
                ROOT.THStack(
                    "hs" + str(channel.GetName()),
                    "hs" + str(channel.GetName()),
                )
            ]
            sample_histograms = []

            original_sample_bin_values = [
                0
            ] * channel.GetData().GetHisto().GetNbinsX()

            for i, sample in enumerate(
                channel.GetSamples()
            ):  # loop over samples to get prefit histograms
                hist = sample.GetHisto()
                hist.SetFillColor(self.predefined_colors[i])
                hist.SetLineColor(ROOT.kBlack)
                sample_histograms += [hist]

                for i in range(hist.GetNbinsX()):
                    original_sample_bin_values[i] += hist.GetBinContent(i + 1)

                self.hs_stacks[-1].Add(
                    sample_histograms[-1]
                )  # save histogram to THStack

            channel_name = "_".join(channel.GetName().split("_")[1:])
            self.hs_stacks[-1].SetTitle(channel_name + " PREFIT")
            self.hs_stacks[-1].Draw("hist")
            self.hs_stacks[-1].GetXaxis().SetLabelSize(0.06)

            maximum_y_val = (
                ROOT.gPad.GetUymax()
            )  # used for moving upper limit of Y ax

            bin_index = 1

            for bin_index in range(1, sample_histograms[-1].GetNbinsX() + 1):
                leftEdge = sample_histograms[-1].GetBinLowEdge(bin_index)
                binWidth = sample_histograms[-1].GetBinWidth(bin_index)
                rightEdge = leftEdge + binWidth

                central_value = original_sample_bin_values[bin_index - 1]
                unc = prefit_unc[bin_index - 1]
                down_value = central_value - unc
                up_value = central_value + unc

                if up_value > maximum_y_val:
                    maximum_y_val = up_value  # if box are over upper limit -> we are going to move upper limit after

                self.boxes += [
                    ROOT.TBox(leftEdge, down_value, rightEdge, up_value)
                ]  # draw uncertainty box for each bin
                self.boxes[-1].SetFillStyle(3004)
                self.boxes[-1].SetFillColor(ROOT.kGray + 3)
                self.boxes[-1].Draw("same")

            self.hs_stacks[-1].SetMaximum(
                1.1 * maximum_y_val
            )  # moving upper limit if required

            self.data_histogram = (
                channel.GetData().GetHisto()
            )  # getting data histogram from channel data
            self.data_histogram.SetStats(0)  # remove extra text over histogram
            self.data_histogram.SetMarkerStyle(3)  # stylization
            self.data_histogram.SetMarkerSize(0.5)

            self.data_plots += [ROOT.TGraphErrors()]
            for i in range(self.data_histogram.GetNbinsX()):
                self.data_plots[-1].SetPoint(
                    i,
                    self.data_histogram.GetBinCenter(i + 1),
                    self.data_histogram.GetBinContent(i + 1),
                )
                self.data_plots[-1].SetPointError(
                    i, 0, self.data_histogram.GetBinError(i + 1) / 2
                )

            self.data_plots[-1].SetMarkerStyle(8)
            self.data_plots[-1].SetMarkerSize(0.5)
            self.data_plots[-1].Draw("same p")  # draw TGraph instead of TH1F

            pad1_bottom.cd()
            number_of_bins = self.data_histogram.GetNbinsX()
            self.bias_graphs += [
                ROOT.TGraphErrors(number_of_bins)
            ]  # TGraph which will be use to save data/bin_value data
            self.bias_graphs[-1].SetTitle("")
            self.bias_graphs[-1].SetMarkerSize(0.4)
            self.bias_graphs[-1].SetMarkerStyle(8)

            for i in range(1, number_of_bins + 1):
                original_value = original_sample_bin_values[i - 1]
                data_value = self.data_histogram.GetBinContent(i)
                data_error = self.data_histogram.GetBinError(i) / 2
                self.bias_graphs[-1].SetPoint(
                    i - 1,
                    self.data_histogram.GetBinCenter(i),
                    data_value / original_value,
                )
                self.bias_graphs[-1].SetPointError(
                    i - 1, 0, data_error / original_value
                )

            self.bias_graphs[-1].Draw("AP")
            self.bias_graphs[-1].GetXaxis().SetLabelSize(0.12)
            self.bias_graphs[-1].GetYaxis().SetLabelSize(0.08)

            for i in range(1, number_of_bins + 1):
                original_value = original_sample_bin_values[i - 1]
                unc = prefit_unc[i - 1]
                up_value = 1 + unc / (original_value)
                down_value = 1 - unc / (original_value)

                if (
                    down_value < 0.5
                ):  # draw limits, by default in ROOT: if down_value << 0.5 -> box will not be drawn
                    down_value = 0.5  # handle it manually
                if (
                    up_value > 1.5
                ):  # by default in ROOT: if up_value >> 1.5 -> box would not be drawn
                    up_value = 1.5

                leftEdge = self.data_histogram.GetBinLowEdge(
                    i
                )  # left edge of bin -> left edge of error box
                binWidth = self.data_histogram.GetBinWidth(i)
                rightEdge = (
                    leftEdge + binWidth
                )  # right edge of bin -> right edge of error box

                self.error_boxes_prefit += [
                    ROOT.TBox(leftEdge, down_value, rightEdge, up_value)
                ]
                self.error_boxes_prefit[-1].SetFillStyle(3004)
                self.error_boxes_prefit[-1].SetFillColor(ROOT.kGray + 3)
                self.error_boxes_prefit[-1].Draw("same")  # drawing error boxes

            minimal_bin_value = self.data_histogram.GetBinLowEdge(1)
            maximum_bin_value = self.data_histogram.GetBinLowEdge(
                number_of_bins
            ) + self.data_histogram.GetBinWidth(number_of_bins)

            self.bias_graphs[-1].GetYaxis().SetRangeUser(0.5, 1.5)
            self.bias_graphs[-1].GetXaxis().SetRangeUser(
                minimal_bin_value, maximum_bin_value
            )

            self.normal_lines += [
                ROOT.TLine(minimal_bin_value, 1, maximum_bin_value, 1)
            ]
            self.normal_lines[-1].SetLineStyle(2)
            self.normal_lines[-1].SetLineWidth(1)
            self.normal_lines[-1].Draw(
                "same"
            )  # normal line, which represent data/bin_value = 1

            pad2_upper.cd()  # same for postfit

            if no_fit:
                (
                    postfit_yields,
                    postfit_yields_uncert,
                    postfit_sample_values,
                ) = self.get_yields_no_fit(
                    obs_var, observables, channel_pdf, result
                )  # if we do for extra workspaces -> use no_fit function
                # which handle all parameters copying
            else:
                (
                    postfit_yields,
                    postfit_yields_uncert,
                    postfit_sample_values,
                ) = self.get_yields(
                    obs_var, observables, channel_pdf, result, prefit=False
                )

            self.second_hs_stacks += [
                ROOT.THStack("fitted stack", "fitted stack")
            ]

            color_number = 0

            for postfit in postfit_sample_values:
                temp_histo = ROOT.TH1F(
                    channel_name + "afterfit" + str(postfit),
                    channel_name + "afterfit" + str(postfit),
                    len(postfit_yields),
                    minimal_bin_value,
                    maximum_bin_value,
                )
                temp_histo.SetFillColor(self.predefined_colors[color_number])
                color_number += 1
                bin_index = 1
                for bin_value in postfit_sample_values[postfit]:
                    temp_histo.SetBinContent(bin_index, bin_value)
                    bin_index += 1
                self.second_hs_stacks[-1].Add(temp_histo)

            channel_name = "_".join(channel.GetName().split("_")[1:])
            self.second_hs_stacks[-1].SetTitle(channel_name + " POSTFIT")
            self.second_hs_stacks[-1].Draw("hist")
            self.second_hs_stacks[-1].GetXaxis().SetLabelSize(0.06)
            self.second_hs_stacks[-1].SetMaximum(1.1 * maximum_y_val)

            bin_index = 1
            for bin_index in range(1, temp_histo.GetNbinsX() + 1):
                leftEdge = temp_histo.GetBinLowEdge(bin_index)
                binWidth = temp_histo.GetBinWidth(bin_index)
                rightEdge = leftEdge + binWidth

                central_value = postfit_yields[bin_index - 1]
                unc = postfit_yields_uncert[bin_index - 1]
                down_value = central_value - unc
                up_value = central_value + unc

                self.boxes += [
                    ROOT.TBox(leftEdge, down_value, rightEdge, up_value)
                ]
                self.boxes[-1].SetFillStyle(3004)
                self.boxes[-1].SetFillColor(ROOT.kGray + 3)
                self.boxes[-1].Draw("same")

            self.data_plots[-1].Draw("same p")

            self.bias_second_graphs += [ROOT.TGraphErrors(number_of_bins)]
            self.bias_second_graphs[-1].SetTitle("")
            self.bias_second_graphs[-1].SetMarkerSize(0.4)
            self.bias_second_graphs[-1].SetMarkerStyle(8)

            pad2_bottom.cd()

            for i in range(1, number_of_bins + 1):
                original_value = postfit_yields[i - 1]
                data_value = self.data_histogram.GetBinContent(i)
                data_error = self.data_histogram.GetBinError(i) / 2

                self.bias_second_graphs[-1].SetPoint(
                    i - 1,
                    self.data_histogram.GetBinCenter(i),
                    data_value / original_value,
                )
                self.bias_second_graphs[-1].SetPointError(
                    i - 1, 0, data_error / original_value
                )

            self.bias_second_graphs[-1].Draw("AP")
            self.bias_second_graphs[-1].GetXaxis().SetLabelSize(0.12)
            self.bias_second_graphs[-1].GetYaxis().SetLabelSize(0.08)

            for i in range(1, number_of_bins + 1):
                original_value = postfit_yields[i - 1]
                unc = postfit_yields_uncert[i - 1]
                up_value = 1 + unc / (original_value)
                down_value = 1 - unc / (original_value)

                leftEdge = self.data_histogram.GetBinLowEdge(i)
                binWidth = self.data_histogram.GetBinWidth(i)
                rightEdge = leftEdge + binWidth

                self.error_boxes += [
                    ROOT.TBox(leftEdge, down_value, rightEdge, up_value)
                ]
                self.error_boxes[-1].SetFillStyle(3004)
                self.error_boxes[-1].SetFillColor(ROOT.kGray + 3)
                self.error_boxes[-1].Draw("same")

            minimal_bin_value = self.data_histogram.GetBinLowEdge(1)
            maximum_bin_value = self.data_histogram.GetBinLowEdge(
                number_of_bins
            ) + self.data_histogram.GetBinWidth(number_of_bins)

            self.bias_second_graphs[-1].GetYaxis().SetRangeUser(0.5, 1.5)
            self.bias_second_graphs[-1].GetXaxis().SetRangeUser(
                minimal_bin_value, maximum_bin_value
            )

            self.normal_lines += [
                ROOT.TLine(minimal_bin_value, 1, maximum_bin_value, 1)
            ]
            self.normal_lines[-1].SetLineStyle(2)
            self.normal_lines[-1].SetLineWidth(1)
            self.normal_lines[-1].Draw("same")

            self.cv_array[-1].SaveAs(channel.GetName() + "_histo.png")
            self.cv_array[-1].Draw()

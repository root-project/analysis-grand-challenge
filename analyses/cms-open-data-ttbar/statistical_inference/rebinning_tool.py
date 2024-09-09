import os

import ROOT


class RebinningTool:
    """
    Rebinning tool, which can be used for rebinning and cutting edges
    Save histograms to new .root file, which have structure as input one

    set_xmin: set lower edge for all histograms
    set_xmax: set upper edge for all histograms
    set_rebin: set rebin coefficient

    example of usage:

    rebinning_tool = RebinningTool()
    rebinning_tool.set_xmin(110) # if xmax was not set -> keep original one
    rebinning_tool.set_rebin(2)
    rebinning_tool.set_input_path("data/histograms.root")
    rebinning_tool.set_output_path("data/temp_histos.root")
    rebinning_tool.apply_rebinning()
    """

    def __init__(self):
        self.fHistoBinLow = None
        self.fHistoBinHigh = None
        self.fRebin = None
        self.fRebinning = False

    def set_xmin(self, xmin):
        self.fHistoBinLow = xmin

    def set_xmax(self, xmax):
        self.fHistoBinHigh = xmax

    def set_rebin(self, rebin):
        self.fRebinning = True
        self.fRebin = rebin

    def list_histograms(self, directory, path=""):
        keys = directory.GetListOfKeys()

        # Loop over all keys
        for key in keys:
            # Get the object associated with the key
            obj = key.ReadObj()

            # Construct the full path of the object

            obj_path = os.path.join(path, obj.GetName())

            # Check if the object is a directory
            if obj.InheritsFrom("TDirectory"):
                # Recursively list histograms in subdirectories
                self.list_histograms(obj, obj_path)
            # Check if the object is a histogram
            elif obj.InheritsFrom("TH1"):
                self.all_histograms.append(obj_path)

    def is_integer(self, value):
        return int(value) == value

    def rebin_histogram(self, original):
        if not self.fRebinning:
            return original

        left_edge = self.fHistoBinLow
        right_edge = self.fHistoBinHigh

        if left_edge is None:
            left_edge = original.GetXaxis().GetXmin()

        if right_edge is None:
            right_edge = original.GetXaxis().GetXmax()

        original_bin_width = original.GetXaxis().GetBinWidth(1)
        assert_check_left = (
            left_edge - original.GetXaxis().GetXmin()
        ) / original_bin_width
        assert_check_right = (
            original.GetXaxis().GetXmax() - right_edge
        ) / original_bin_width

        if not self.is_integer(assert_check_left) or not self.is_integer(assert_check_right):
            print("Error: The left_edge and right_edge are not multiples of the original bin width")
            return original

        number_of_remove_bins = (
            (left_edge - original.GetXaxis().GetXmin())
            + (original.GetXaxis().GetXmax() - right_edge)
        ) / original_bin_width
        new_nbins = int(original.GetNbinsX() - number_of_remove_bins)

        original_new = ROOT.TH1F(
            original.GetName() + "temp_rebin_clone",
            original.GetTitle(),
            new_nbins,
            left_edge,
            right_edge,
        )
        skipped_bins_left = int((left_edge - original.GetXaxis().GetXmin()) / original_bin_width)

        for i in range(1, new_nbins + 1):
            bin_idx = i + skipped_bins_left
            original_new.SetBinContent(i, original.GetBinContent(bin_idx))
            original_new.SetBinError(i, original.GetBinError(bin_idx))

        output = original_new.Rebin(self.fRebin)
        output.SetDirectory(ROOT.nullptr)

        return output

    def apply_rebinning(self, input_file_path=None, output_file_path=None):
        if input_file_path is None:
            input_file_path = self.input_path
        if output_file_path is None:
            output_file_path = self.output_path

        file = ROOT.TFile(input_file_path, "READ")

        output_file = ROOT.TFile(output_file_path, "RECREATE")
        output_file.Close()

        self.all_histograms = []
        self.list_histograms(file)

        for hist_path in self.all_histograms:
            hist = file.Get(hist_path)
            hist_rebinned = self.rebin_histogram(hist)
            output_file = ROOT.TFile(output_file_path, "UPDATE")

            if "/" in hist_path:
                dir_path = hist_path.rsplit("/", 1)[0]
                self.mkdir_p(output_file, dir_path)
                output_file.cd(dir_path)

            hist_rebinned.Write(hist_path.split("/")[-1])
            output_file.Close()

        file.Close()

    def set_input_path(self, input_path):
        self.input_path = input_path

    def set_output_path(self, output_path):
        self.output_path = output_path
        file = ROOT.TFile(self.output_path, "RECREATE")
        file.Close()

import ROOT


class AGCSample():
    """
    class created to cover extra functions, which are implemented in cabinetry, but missed in cabinetry
    """

    def __init__(self, name, histo_name, histo_file, histo_path=""):
        self.root_sample = ROOT.RooStats.HistFactory.Sample(name, histo_name, histo_file, histo_path)
        self.output_path = "statistical_data/HistFactoryExtra.root"  
        # since ROOT need to collect histogram from file at some point, we store created histograms in additional file
        # probably can be changed in next ROOT releases -> histograms can be stored just in RAM
        self.fInputFile = ""  # default input file for all systematics

    # set same input file for all provided systematics: HistoSysm, HistoFactor, NormPlusShape
    def SetSystematicsInputFile(self, file):  
        self.fInputFile = file

    def GetHistFactorySample(self):
        return self.root_sample
    
    def AddHistoSys(
        self,
        name,
        histoname_up=None,
        histofile_up=None,
        histopath_up="",  # overrided function, allowing not to provide additional input file
        histoname_down=None,
        histofile_down=None,
        histopath_down="",
    ):

        histSys = ROOT.RooStats.HistFactory.HistoSys()

        if histofile_up is None:
            histofile_down = self.fInputFile
            histofile_up = self.fInputFile

        histSys.SetName(name)
        histSys.SetHistoName(self.GetHistFactorySample().GetHistoName())
        histSys.SetHistoPathHigh(histopath_up)
        histSys.SetHistoPathLow(histopath_down)
        histSys.SetHistoNameHigh(histoname_up)
        histSys.SetHistoNameLow(histoname_down)
        histSys.SetInputFileHigh(histofile_up)
        histSys.SetInputFileLow(histofile_down)

        self.GetHistFactorySample().GetHistoSysList().push_back(histSys)

    def AddHistoFactor(
        self,
        name,
        histoname_up=None,
        histofile_up=None,
        histopath_up="",  # overrided function, allowing not to provide additional input file
        histoname_down=None,
        histofile_down=None,
        histopath_down="",
    ):

        histFactor = ROOT.RooStats.HistFactory.HistoFactor()

        if histofile_up is None:
            histofile_down = self.fInputFile
            histofile_up = self.fInputFile

        histFactor.SetName(name)
        histFactor.SetHistoName(self.GetHistFactorySample().GetHistoName())
        histFactor.SetHistoPathHigh(histopath_up)
        histFactor.SetHistoPathLow(histopath_down)
        histFactor.SetHistoNameHigh(histoname_up)
        histFactor.SetHistoNameLow(histoname_down)
        histFactor.SetInputFileHigh(histofile_up)
        histFactor.SetInputFileLow(histofile_down)

        self.GetHistFactorySample().GetHistoFactorList().push_back(histFactor)

    def AddShapeSys(
        self,
        name,
        constraint_type=None,
        histoname=None,
        histofile=None,
        histopath="",
    ):  # overrided function, allowing not to provide additional input file

        shapeSys = ROOT.RooStats.HistFactory.ShapeSys()

        if histofile is None:
            histofile = self.fInputFile

        shapeSys.SetName(name)
        shapeSys.SetHistoName(self.GetHistFactorySample().GetHistoName())
        shapeSys.SetHistoPath(histopath)
        shapeSys.SetInputFile(histofile)

        shapeSys.SetConstraintType(constraint_type)

        self.GetHistFactorySample().GetShapeSysList().push_back(shapeSys)

    def AddNormPlusShapeHistoSys(
        self,
        name,
        histoname_up=None,
        histofile_up=None,
        histopath_up="",  # check more here: https://github.com/scikit-hep/cabinetry/issues/26
        histoname_down=None,
        histofile_down=None,
        histopath_down="",
    ):
        if histofile_up is None:
            histofile_up = self.fInputFile
            histofile_down = self.fInputFile
            assert (histofile_up != ""), "ERROR: You not specified input file for sample"

        if histoname_down is None:
            self.Symmetrize_AddNormPlusShapeHistoSys(name, histoname_up, histofile_up, histopath_up)
        else:
            self.NonSymmetrize_AddNormPlusShapeHistoSys(
                name,
                histoname_up,
                histofile_up,
                histopath_up,
                histoname_down,
                histofile_down,
                histopath_down,
            )

    def Symmetrize_AddNormPlusShapeHistoSys(
        self, name, histoname, histofile, histopath
    ):
        """
        first symmertrize histogram -> calculate overallsys -> save all to sample
        check more here:
        https://github.com/scikit-hep/cabinetry/issues/26
        """
        channel_name = str(histoname.split("_")[0])

        file = ROOT.TFile(histofile, "READ")
        dir = file.GetDirectory(histopath)
        hist_top = dir.Get(histoname)

        hist_nominal_file = ROOT.TFile(self.GetHistFactorySample().GetInputFile(), "READ")
        hist_nominal_name = self.GetHistFactorySample().GetHistoName()
        hist_nominal_directory = hist_nominal_file.GetDirectory(self.GetHistFactorySample().GetHistoPath())
        hist_nominal = hist_nominal_directory.Get(hist_nominal_name)

        norm_factor_up = hist_top.Integral() / hist_nominal.Integral()
        h_new = hist_top.Clone(f"{channel_name}_{self.GetHistFactorySample().GetName()}_{name}_norm_plus_shape_up_clone")
        h_new.Scale(1 / norm_factor_up)

        h_down = hist_nominal.Clone(f"{channel_name}_{self.GetHistFactorySample().GetName()}_{name}_norm_plus_shape_down_clone")
        h_down.Scale(2)
        h_down.Add(h_new, -1)

        output_file = ROOT.TFile(self.output_path, "UPDATE")

        hist_up_name = str(f"{channel_name}_{self.GetHistFactorySample().GetName()}_{name}_norm_plus_shape_up")
        hist_down_name = str(f"{channel_name}_{self.GetHistFactorySample().GetName()}_{name}_norm_plus_shape_down")

        h_new.Write(hist_up_name)
        h_down.Write(hist_down_name)

        output_file.Close()

        histSys = ROOT.RooStats.HistFactory.HistoSys()

        histSys.SetName(name)
        histSys.SetHistoPathHigh("")
        histSys.SetHistoPathLow("")

        histSys.SetHistoNameHigh(hist_up_name)
        histSys.SetHistoNameLow(hist_down_name)

        histSys.SetInputFileLow(self.output_path)
        histSys.SetInputFileHigh(self.output_path)

        overallSys = ROOT.RooStats.HistFactory.OverallSys()
        overallSys.SetName(name)
        overallSys.SetLow(2 - norm_factor_up)
        overallSys.SetHigh(norm_factor_up)

        self.GetHistFactorySample().GetHistoSysList().push_back(histSys)
        self.GetHistFactorySample().GetOverallSysList().push_back(overallSys)

    def NonSymmetrize_AddNormPlusShapeHistoSys(
        self,
        name,
        histoname_up,
        histofile_up,
        histopath_up,
        histoname_down,
        histofile_down,
        histopath_down,
    ):
        channel_name = str(histoname_up.split("_")[0])

        file = ROOT.TFile(histofile_up, "READ")
        dir = file.GetDirectory(histopath_up)
        hist_top = dir.Get(histoname_up)

        hist_nominal_file = ROOT.TFile(self.GetHistFactorySample().GetInputFile(), "READ")
        hist_nominal_name = self.GetHistFactorySample().GetHistoName()
        hist_nominal_directory = hist_nominal_file.GetDirectory(self.GetHistFactorySample().GetHistoPath())
        hist_nominal = hist_nominal_directory.Get(hist_nominal_name)

        norm_factor_up = hist_top.Integral() / hist_nominal.Integral()
        h_new_up = hist_top.Clone(f"{channel_name}_{self.GetHistFactorySample().GetName()}_{name}_norm_plus_shape_up_clone")
        h_new_up.Scale(1 / norm_factor_up)

        file_down = ROOT.TFile(histofile_down, "READ")
        dir_down = file_down.GetDirectory(histopath_down)
        hist_down = dir_down.Get(histoname_down)

        norm_factor_down = hist_down.Integral() / hist_nominal.Integral()
        h_new_down = hist_down.Clone(f"{channel_name}_{self.GetHistFactorySample().GetName()}_{name}_norm_plus_shape_down_clone")
        h_new_down.Scale(1 / norm_factor_down)

        output_file = ROOT.TFile(self.output_path, "UPDATE")

        hist_up_name = str(f"{channel_name}_{self.GetHistFactorySample().GetName()}_{name}_norm_plus_shape_up")
        hist_down_name = str(f"{channel_name}_{self.GetHistFactorySample().GetName()}_{name}_norm_plus_shape_down")

        h_new_up.Write(hist_up_name)
        h_new_down.Write(hist_down_name)

        output_file.Close()

        histSys = ROOT.RooStats.HistFactory.HistoSys()

        histSys.SetName(name)
        histSys.SetHistoPathHigh("")
        histSys.SetHistoPathLow("")

        histSys.SetHistoNameHigh(hist_up_name)
        histSys.SetHistoNameLow(hist_down_name)

        histSys.SetInputFileLow(self.output_path)
        histSys.SetInputFileHigh(self.output_path)

        overallSys = ROOT.RooStats.HistFactory.OverallSys()
        overallSys.SetName(name)
        overallSys.SetLow(norm_factor_down)
        overallSys.SetHigh(norm_factor_up)

        self.GetHistFactorySample().GetHistoSysList().push_back(histSys)
        self.GetHistFactorySample().GetOverallSysList().push_back(overallSys)
        
    def ActivateStatError(self):
        self.GetHistFactorySample().ActivateStatError()
        
    def AddOverallSys(self, overallsys):
        self.GetHistFactorySample().AddOverallSys(overallsys)
        
    def AddNormFactor(self, name, central, left, right):
        self.GetHistFactorySample().AddNormFactor(name, central, left, right)

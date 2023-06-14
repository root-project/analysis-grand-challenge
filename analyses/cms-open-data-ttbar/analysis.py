import argparse
import json
import os
import time

from urllib.request import urlretrieve

PARSER = argparse.ArgumentParser()
PARSER.add_argument("--ncores",
                    "-c",
                    help=("How many cores to use. If choosing a distributed execution, "
                          "this is the amount of cores per node."),
                    default = len(os.sched_getaffinity(0)),
                    type=int)
PARSER.add_argument("--scheduling-mode",
                    "-s",
                    help=("The scheduling mode of the analysis. RDataFrame supports "
                          "both single-node and multi-node parallelization."),
                    default="imt",
                    choices=["imt", "dask-local", "dask-ssh"])
PARSER.add_argument("--nodes", help="String containing the list of hostnames to be used. Useful only in pair with 'dask-ssh'"
                                    " choice for the '--scheduling-mode' parameter",
                    type=str, default=None)
PARSER.add_argument("--npartitions", help="How many partitions to use.", type=int, default=None)
PARSER.add_argument("--n-files-max-per-sample", "-f",
                    help="How many files per sample will be processed. Default -1 (all files for all samples).",
                    type=int, default=-1)
PARSER.add_argument("--storage-location", "-l",
                    help="Where the data resides. Default downloads original dataset stored at UNL.",
                    default="unl", choices=["unl", "cern-xrootd"])
PARSER.add_argument("--histograms-output-file",
                    help="Name of the output file to store histograms.", default="histograms.root")
PARSER.add_argument("--download", "-d", help="Download the files locally when executing.", action="store_true")
PARSER.add_argument("-v", "--verbose", action="store_true")
ARGS = PARSER.parse_args()

# importing ROOT takes a little while, so we only do it if argument parsing succeeded
import ROOT # noqa: {E402}

if ARGS.scheduling_mode == "imt":

    RDataFrame = ROOT.RDataFrame
    RunGraphs = ROOT.RDF.RunGraphs
    VariationsFor = ROOT.RDF.Experimental.VariationsFor

    def init_functions():
        ROOT.gSystem.CompileMacro("helper.cpp", "kO")

else:

    # Dask configuration useful in distributed mode
    RDataFrame = ROOT.RDF.Experimental.Distributed.Dask.RDataFrame
    RunGraphs = ROOT.RDF.Experimental.Distributed.RunGraphs
    VariationsFor = ROOT.RDF.Experimental.Distributed.VariationsFor
    initialize = ROOT.RDF.Experimental.Distributed.initialize

    if ARGS.scheduling_mode == "dask-ssh":
        if not ARGS.nodes:
            raise ValueError("For SSHCluster deployments, please specify a "
                             "string with a comma-separated list of hostnames of nodes that will be used.")

        if ARGS.npartitions is None:
            n_compute_nodes = len(ARGS.nodes.split(",")) - 1
            ARGS.npartitions = ARGS.ncores * n_compute_nodes

    elif ARGS.scheduling_mode == "dask-local":
        if ARGS.npartitions is None:
            ARGS.npartitions = ARGS.ncores

    from distributed import Client, LocalCluster, SSHCluster, get_worker

    def create_localcluster_connection(ncores: int) -> Client:
        cluster = LocalCluster(n_workers=ncores, threads_per_worker=1, processes=True)
        client = Client(cluster)
        return client

    def create_sshcluster_connection(nodes: str, ncores: int) -> Client:
        parsed_nodes = nodes.split(',')
        scheduler = parsed_nodes[:1]
        workers = parsed_nodes[1:]

        print(f"List of nodes: {scheduler=}, {workers=}")

        # The creation of the SSHCluster object needs to be further configured according to needs.
        # For example, in some clusters the "local_directory" key must be supplied in the worker_options dictionary.
        cluster = SSHCluster(scheduler + workers,
                             connect_options={"known_hosts": None},
                             worker_options={"nprocs": ncores, "nthreads": 1, "memory_limit": "32GB"})

        return Client(cluster)

    def create_connection(nodes: str, ncores: int, scheduling_mode: str) -> Client:
        if scheduling_mode == "dask-local":
            return create_localcluster_connection(ncores)
        elif scheduling_mode == "dask-ssh":
            return create_sshcluster_connection(nodes, ncores)
        # Add more cluster types here to accomodate different deployments
        else:
            raise ValueError(
                f"Unexpected scheduling mode '{scheduling_mode}'. Acceptable values are ['dask-local', 'dask-ssh'].")

    def init_functions():
        try:
            localdir = get_worker().local_directory
            helper_path = os.path.join(localdir, "helper.cpp")
        except ValueError:
            # get_worker raises an error in case it is called from the local machine
            # for now work around this by silencing the error.
            helper_path = "helper.cpp"

        ROOT.gSystem.CompileMacro(helper_path, "kO")

if ARGS.verbose:
    verbosity = ROOT.Experimental.RLogScopedVerbosity(
        ROOT.Detail.RDF.RDFLogChannel(), ROOT.Experimental.ELogLevel.kInfo)


class TtbarAnalysis(dict):

    def __init__(self, n_files_max_per_sample, download_input_data, storage_location, num_bins=25, bin_low=50, bin_high=550, connection=None):

        # Store input arguments
        self.n_files_max_per_sample = n_files_max_per_sample  # the number of files to be processed per sample
        self.download_input_data = download_input_data
        self.storage_location = storage_location
        self.ntuples_file = "ntuples.json"
        self.num_bins = num_bins
        self.bin_low = bin_low
        self.bin_high = bin_high

        # Connection handle in case distributed mode was selected
        self.connection = connection

        self.variations = {}  # serves as temporary storage for all histograms produced by VariationsFor
        self._nevts_total = {}
        # dictionary assigning file URLs (paths) to each process, variation, and region
        self.input_data = self._construct_fileset()
        # using https://atlas-groupdata.web.cern.ch/atlas-groupdata/dev/AnalysisTop/TopDataPreparation/XSection-MC15-13TeV.data
        # for reference
        # x-secs are in pb
        self.xsec_info = {
            "ttbar": 396.87 + 332.97,  # nonallhad + allhad, keep same x-sec for all
            "single_top_s_chan": 2.0268 + 1.2676,
            "single_top_t_chan": (36.993 + 22.175) / 0.252,  # scale from lepton filter to inclusive
            "single_top_tW": 37.936 + 37.906,
            "wjets": 61457 * 0.252,  # e/mu+nu final states
            "data": None
        }

    def _optionally_download_data(self, file_paths, process, variation):
        if (self.download_input_data):
            dir_name = f"input/{process}_{variation}"
            os.makedirs(dir_name, exist_ok=True)
            for i in range(len(file_paths)):
                path = file_paths[i]
                file = f"{dir_name}/{i}.root"
                if not os.path.exists(file):
                    urlretrieve(path, file)
                    print(f"{file} has been created")
                else:
                    print(f"{file} already exists")

    def _construct_fileset(self):

        with open(self.ntuples_file) as f:
            file_info = json.load(f)

        fileset = {}
        for process in file_info.keys():
            if process == "data":
                continue  # skip data
            fileset[process] = {}
            self[process] = {}
            self._nevts_total[process] = {}

            for variation in file_info[process].keys():
                file_list = file_info[process][variation]["files"]
                if self.n_files_max_per_sample != -1:
                    file_list = file_list[:self.n_files_max_per_sample]  # use partial set of samples
                file_paths = [f["path"] for f in file_list]
                if (self.storage_location == "cern-xrootd"):
                    file_paths = [f.replace("https://xrootd-local.unl.edu:1094//store/user/AGC",
                                            "root://eoscms.cern.ch//eos/cms/store/test/agc") for f in file_paths]

                fileset[process].update({variation: file_paths})
                nevts_total = sum([f["nevts"] for f in file_list])
                self._nevts_total[process].update({variation: nevts_total})
                self[process][variation] = {}

                self._optionally_download_data(file_paths, process, variation)

        return fileset

    def fill(self, process: str, variation: str):

        # all operations are handled by RDataFrame class, so the first step is the RDataFrame object instantiating
        input_data = self.input_data[process][variation]
        if ARGS.scheduling_mode == "imt":
            d = RDataFrame("events", input_data)
        else:
            d = RDataFrame("events", input_data, daskclient=self.connection, npartitions=ARGS.npartitions)
            d._headnode.backend.distribute_unique_paths(["helper.cpp", ])

        # normalization for MC
        x_sec = self.xsec_info[process]
        nevts_total = self._nevts_total[process][variation]
        lumi = 3378  # /pb
        xsec_weight = x_sec * lumi / nevts_total
        d = d.Define("weights", str(xsec_weight))  # default weights

        if variation == "nominal":

            # jet_pt variations definition
            # pt_scale_up() and pt_res_up(jet_pt) return scaling factors applying to jet_pt
            # pt_scale_up() - jet energy scaly systematic
            # pt_res_up(jet_pt) - jet resolution systematic

            d = d.Vary("jet_pt",
                       "ROOT::RVec<ROOT::RVecF>{jet_pt*pt_scale_up(), jet_pt*jet_pt_resolution(jet_pt.size())}",
                       ["pt_scale_up", "pt_res_up"])
            if process == "wjets":

                # flat weight variation definition
                d = d.Vary("weights",
                           "weights*flat_variation()",
                           [f"scale_var_{direction}" for direction in ["up", "down"]]
                           )

        # event selection - the core part of the algorithm applied for both regions
        # selecting events containing at least one lepton and four jets with pT > 25 GeV
        # applying requirement at least one of them must be b-tagged jet (see details in the specification)
        d = d.Define("electron_pt_mask", "electron_pt>25").Define("muon_pt_mask", "muon_pt>25").Define("jet_pt_mask", "jet_pt>25")\
             .Filter("Sum(electron_pt_mask) + Sum(muon_pt_mask) == 1")\
             .Filter("Sum(jet_pt_mask) >= 4")\
             .Filter("Sum(jet_btag[jet_pt_mask]>=0.5)>=1")

        # b-tagging variations for nominal samples
        d = d.Vary("weights",
                   "ROOT::RVecD{weights*btag_weight_variation(jet_pt[jet_pt_mask])}",
                   [f"{weight_name}_{direction}" for weight_name in [f"btag_var_{i}" for i in range(4)] for direction in [
                       "up", "down"]]
                   ) if variation == "nominal" else d

        # as next steps for different regions are different, there is a fork in the algorithm
        # we create another RDF pointer for each region called "fork"
        measured = {"4j1b": "HT", "4j2b": "trijet_mass"}  # columns names of observables for two regions
        for region in ["4j1b", "4j2b"]:
            observable = measured[region]

            if region == "4j1b":

                # only one b-tagged region required
                # observable is total transvesre momentum
                fork = d.Filter("Sum(jet_btag[jet_pt_mask]>=0.5)==1").Define(observable, "Sum(jet_pt[jet_pt_mask])")

            elif region == "4j2b":

                # select events with at least 2 b-tagged jets
                # building four-momentum vectors for each jet
                fork = (
                    d.Filter("Sum(jet_btag[jet_pt_mask]>=0.5)>1")
                    .Define("jet_p4",
                            """
                            ROOT::VecOps::Construct<ROOT::Math::PxPyPzMVector>(
                                jet_px[jet_pt_mask], jet_py[jet_pt_mask], jet_pz[jet_pt_mask], jet_mass[jet_pt_mask])""")
                )
                # building trijet combinations
                fork = fork.Define("trijet",
                                   "ROOT::VecOps::Combinations(jet_pt[jet_pt_mask],3)"
                                   ).Define("ntrijet", "trijet[0].size()")

                # assigning four-momentums to each trijet combination
                fork = fork.Define("trijet_p4",
                                   """
                                   ROOT::RVec<ROOT::Math::PxPyPzMVector> trijet_p4(ntrijet);
                                   for (int i = 0; i < ntrijet; ++i)
                                   {
                                       int j1 = trijet[0][i];
                                       int j2 = trijet[1][i];
                                       int j3 = trijet[2][i];
                                       trijet_p4[i] = jet_p4[j1] + jet_p4[j2] + jet_p4[j3];
                                   }
                                   return trijet_p4;
                                   """
                                   )

                # getting trijet transverse momentum values from four-momentum vectors
                fork = fork.Define("trijet_pt",
                                   "return ROOT::VecOps::Map(trijet_p4, [](ROOT::Math::PxPyPzMVector v) { return v.Pt(); })"
                                   )

                # trijet_btag is a helpful array of bool values indicating whether or not the maximum btag value in trijet is larger than 0.5 threshold
                fork = fork.Define("trijet_btag",
                                   """
                                   ROOT::RVecB btag(ntrijet);
                                   for (int i = 0; i < ntrijet; ++i)
                                   {
                                       int j1 = trijet[0][i];
                                       int j2 = trijet[1][i];
                                       int j3 = trijet[2][i];
                                       btag[i] = std::max({jet_btag[j1], jet_btag[j2], jet_btag[j3]}) > 0.5;
                                   }
                                   return btag;
                                   """
                                   )
                # find trijet with maximum pt and higher that threshold btag
                # get mass for found jet four-vector
                # trijet mass themself is an observable quantity
                fork = fork.Define(observable,
                                   """
                                   double mass{};
                                   double Pt{};
                                   double indx{};
                                   for (int i = 0; i < ntrijet; ++i) {
                                       if ((Pt < trijet_pt[i]) && (trijet_btag[i])) {
                                           Pt = trijet_pt[i];
                                           indx = i;
                                       }
                                   }
                                   mass = trijet_p4[indx].M();
                                   return mass;
                                   """
                                   )

            # fill histogram for observable column in RDF object
            res = fork.Histo1D((f"{process}_{variation}_{region}", process, self.num_bins,
                               self.bin_low, self.bin_high), observable, "weights")
            self.hist.append(res)  # save the pointer to further triggering
            print(f"histogram {process}_{variation}_{region} has been created")

            # save pointers for variations
            # self.variations is a temporary container for all pointers
            if variation == "nominal":
                self.variations[f"{process}__{region}"] = VariationsFor(res)
            else:
                self[process][variation][region] = res

    # build 9 Graphs for each data sample
    def Fill(self):
        self.hist = []
        for process in self:
            for variation in self.input_data[process]:
                self.fill(process, variation)

    # run 9 Graphs for each data sample
    def Accumulate(self):
        RunGraphs(self.hist)

    # transform TtbarAnalysis to dictionary (process, variation, region) -> histogram
    def TransfToDict(self):
        for key in self.variations.keys():
            hist_map = self.variations[key]
            key = str(key).split("__")
            process = key[0]
            region = key[1]
            for hist_name in hist_map.GetKeys():
                variation = "nominal" if hist_name == "nominal" else str(hist_name).split(":")[1]
                if variation not in self[process]:
                    self[process][variation] = {}
                hist = hist_map[hist_name]
                if not isinstance(hist, ROOT.TH1D):
                    hist = hist.GetValue()
                self[process][variation][region] = hist
        self.ExportJSON()

    def GetProcStack(self, region, variation="nominal"):
        return [self[process][variation][region] for process in self]

    def GetVarStack(self, region, process="ttbar", variations=None):
        variations = variations if variations is not None else self[process].keys()
        histos = [self[process][variation][region] for variation in variations]
        return [h.GetValue() if not isinstance(h, ROOT.TH1D) else h for h in histos]

    # necessary only for sanity checks
    def ExportJSON(self):
        data = {}
        for process in self:
            data[process] = {}
            for variation in self[process]:
                data[process][variation] = [region for region in self[process][variation]]
        with open("data.json", "w") as f:
            json.dump(data, f)


def analyse(connection=None):

    analysisManager = TtbarAnalysis(download_input_data=ARGS.download,
                                    n_files_max_per_sample=ARGS.n_files_max_per_sample,
                                    storage_location=ARGS.storage_location,
                                    connection=connection)

    # At this stage, analysisManager keeps all file URLs:
    print(f"processes in fileset: {list(analysisManager.keys())}")
    print(
        f'\nexample of information inside analysisManager:\n{{\n  "urls": [{analysisManager.input_data["ttbar"]["nominal"][0]}, ...],')

    t0 = time.time()
    analysisManager.Fill()
    t1 = time.time()

    print(f"\npreprocessing took {round(t1 - t0,2)} seconds")

    analysisManager.Accumulate()
    t2 = time.time()

    print(f"processing took {round(t2 - t1,2)} seconds")
    print(f"execution took {round(t2 - t0,2)} seconds")

    analysisManager.TransfToDict()

    return analysisManager


def make_plots(analysisManager):

    width = 2160
    height = 2160
    c = ROOT.TCanvas("c", "c", width, height)
    ROOT.gStyle.SetPalette(ROOT.kRainBow)

    # Region 1 stack
    hlist = analysisManager.GetProcStack(region="4j1b")
    hs = ROOT.THStack("j4b1", ">=4 jets, 1 b-tag; H_{T} [GeV]")
    for h in hlist:
        h = ROOT.Slice(h, 120, 550)
        ptr = h.Rebin(2, h.GetTitle())
        hs.Add(ptr)
    hs.Draw("hist pfc plc")
    c.Draw()
    x = hs.GetXaxis()
    x.SetTitleOffset(1.5)
    x.CenterTitle()
    c.BuildLegend(0.65, 0.7, 0.9, 0.9)
    c.SaveAs("reg1.png")

    # Region 2 stack
    hlist = analysisManager.GetProcStack(region="4j2b")
    hs = ROOT.THStack("j4b1", ">=4 jets, 2 b-tag; H_{T} [GeV]")
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
    btag_variations = ["nominal", "btag_var_0_up", "btag_var_1_up", "btag_var_2_up", "btag_var_3_up"]
    freshstack = analysisManager.GetVarStack(region="4j1b", variations=btag_variations)

    hs = ROOT.THStack("j4b1btag", "btag-variations ; H_{T} [GeV]")
    for h, name in zip(freshstack, btag_variations):
        print(name)
        ptr = h.Rebin(2, name)
        ptr.SetLineWidth(4)
        ptr.SetTitle(name)
        hs.Add(ptr)
    hs.Draw("hist nostack plc")
    c.Draw()
    x = hs.GetXaxis()
    x.SetRangeUser(120, 500)
    x.SetTitleOffset(1.5)
    x.CenterTitle()
    c.BuildLegend(0.65, 0.7, 0.9, 0.9)
    c.SaveAs("btag.png")

    # Jet energy variations
    jet_variations = ["nominal", "pt_scale_up", "pt_res_up"]
    freshstack = analysisManager.GetVarStack(region="4j2b", variations=jet_variations)
    hs = ROOT.THStack("4j2bjet", "Jet energy variations ; m_{bjj} [GeV]")
    for h, name in zip(freshstack, jet_variations):
        print(name)
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

    # Save histograms to disk
    with ROOT.TFile.Open(ARGS.histograms_output_file, "RECREATE") as output:
        for process in analysisManager:
            for variation in analysisManager[process]:
                for region in analysisManager[process][variation]:
                    hist_name = (
                        f"{region}_{process}_{variation}" if variation != "nominal" else f"{region}_{process}"
                    )
                    hist = analysisManager[process][variation][region]
                    if not isinstance(hist, ROOT.TH1D):
                        hist = hist.GetValue()
                    if hist.IsZombie():
                        raise TypeError(hist_name)
                    hist_sliced = ROOT.Slice(hist, 120, 550)
                    hist_binned = hist_sliced.Rebin(2, hist.GetTitle())
                    output.WriteObject(hist_binned, hist_name)


def main():

    if ARGS.scheduling_mode == "imt":
        init_functions()
        ROOT.EnableImplicitMT(ARGS.ncores)
        print(f"The num of threads = {ROOT.GetThreadPoolSize()}")
        # No handle needed in local mode
        connection = None
    else:
        initialize(init_functions)
        # Create connection to the cluster in distributed mode
        connection = create_connection(ARGS.nodes, ARGS.ncores, ARGS.scheduling_mode)

    results = analyse(connection)
    make_plots(results)

    if connection is not None:
        connection.close()


if __name__ == "__main__":
    raise SystemExit(main())

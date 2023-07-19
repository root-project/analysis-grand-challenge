import argparse
import os
from pathlib import Path
from time import time
from typing import Optional

import ROOT
from distributed import Client, LocalCluster, SSHCluster, get_worker
from plotting import save_plots
from utils import (
    AGCInput,
    AGCResult,
    postprocess_results,
    retrieve_inputs,
    save_histos,
)

# Using https://atlas-groupdata.web.cern.ch/atlas-groupdata/dev/AnalysisTop/TopDataPreparation/XSection-MC15-13TeV.data
# as a reference. Values are in pb.
XSEC_INFO = {
    "ttbar": 396.87 + 332.97,  # nonallhad + allhad, keep same x-sec for all
    "single_top_s_chan": 2.0268 + 1.2676,
    "single_top_t_chan": (36.993 + 22.175) / 0.252,  # scale from lepton filter to inclusive
    "single_top_tW": 37.936 + 37.906,
    "wjets": 61457 * 0.252,  # e/mu+nu final states
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument(
        "--n-max-files-per-sample",
        "-n",
        help="How many files per sample will be processed (if absent, all files for all samples).",
        type=int,
    )
    p.add_argument(
        "--data-cache",
        "-d",
        help="Use the specified directory as a local data cache: required input datasets will be downloaded here and the analysis will read this local copy of the data.",
    )
    p.add_argument(
        "--remote-data-prefix",
        help="""The original data is stored at 'https://xrootd-local.unl.edu:1094//store/user/AGC'.
                If this option is passed, that prefix is replaced with the argument to this option when accessing
                remote data. For example for the version of the input datasets stored on EOS use
                `--remote-data-prefix='root://eoscms.cern.ch//eos/cms/store/test/agc'`.""",
    )
    p.add_argument(
        "--output",
        "-o",
        help="Name of the file where analysis results will be stored. If it already exists, contents are overwritten.",
        default="histograms.root",
    )
    p.add_argument(
        "--scheduler",
        "-s",
        help="""The scheduler for RDataFrame parallelization. Will honor the --ncores flag.
                The default is `mt`, i.e. multi-thread execution.
                If dask-ssh, a list of worker node hostnames to connect to should be provided via the --nodes option.""",
        default="mt",
        choices=["mt", "dask-local", "dask-ssh"],
    )
    p.add_argument(
        "--ncores",
        "-c",
        help=(
            "Number of cores to use. In case of distributed execution this is the amount of cores per node."
        ),
        default=len(os.sched_getaffinity(0)),
        type=int,
    )
    p.add_argument(
        "--npartitions",
        help="Number of data partitions. Only used in case of distributed execution. By default it follows RDataFrame defaults.",
        type=int,
    )
    p.add_argument(
        "--hosts",
        help="A comma-separated list of worker node hostnames. Only required if --scheduler=dask-ssh, ignored otherwise.",
    )
    p.add_argument("-v", "--verbose", help="Turn on verbose execution logs.", action="store_true")

    return p.parse_args()


def create_dask_client(scheduler: str, ncores: int, hosts: str) -> Client:
    """Create a Dask distributed client."""
    if scheduler == "dask-local":
        lc = LocalCluster(n_workers=ncores, threads_per_worker=1, processes=True)
        return Client(lc)

    if scheduler == "dask-ssh":
        workers = hosts.split(",")
        print(f"Using worker nodes: {workers=}")
        # The creation of the SSHCluster object might need to be further configured to fit specific use cases.
        # For example, in some clusters the "local_directory" key must be supplied in the worker_options dictionary.
        sshc = SSHCluster(
            workers,
            connect_options={"known_hosts": None},
            worker_options={"nprocs": ncores, "nthreads": 1, "memory_limit": "32GB"},
        )
        return Client(sshc)

    raise ValueError(
        f"Unexpected scheduling mode '{scheduler}'. Valid modes are ['dask-local', 'dask-ssh']."
    )


def make_rdf(
    files: list[str], client: Optional[Client], npartitions: Optional[int]
) -> ROOT.RDataFrame:
    """Construct and return a dataframe or, if a dask client is present, a distributed dataframe."""
    if client is not None:
        d = ROOT.RDF.Experimental.Distributed.Dask.RDataFrame(
            "Events", files, daskclient=client, npartitions=npartitions
        )
        d._headnode.backend.distribute_unique_paths(
            [
                "helpers.cpp",
            ]
        )
        return d

    return ROOT.RDataFrame("Events", files)


def define_trijet_mass(df: ROOT.RDataFrame) -> ROOT.RDataFrame:
    """Add the trijet_mass observable to the dataframe after applying the appropriate selections."""

    # First, select events with at least 2 b-tagged jets
    df = df.Filter("Sum(Jet_btagCSVV2[Jet_pt_mask]>0.5)>1")

    # Build four-momentum vectors for each jet
    df = (  
        df.Define(
        "Jet_p4",
        """
        ROOT::VecOps::Construct<ROOT::Math::PxPyPzMVector>(
            ROOT::VecOps::Construct<ROOT::Math::PtEtaPhiMVector>(
                Jet_pt[Jet_pt_mask], Jet_eta[Jet_pt_mask], Jet_phi[Jet_pt_mask], Jet_mass[Jet_pt_mask]
            )
        )
        """,
        )
    )

    # Build trijet combinations
    df = df.Define("Trijet", "ROOT::VecOps::Combinations(Jet_pt[Jet_pt_mask],3)")

    
    # Trijet_btag is a helpful array mask indicating whether or not the maximum btag value in Trijet is larger than the 0.5 threshold
    df = df.Define(
            "Trijet_btag",
            """
            auto Jet_btagCSVV2_masked = Jet_btagCSVV2[Jet_pt_mask];
            auto J1_btagCSVV2 = ROOT::VecOps::Take(Jet_btagCSVV2_masked, Trijet[0]);
            auto J2_btagCSVV2 = ROOT::VecOps::Take(Jet_btagCSVV2_masked, Trijet[1]);
            auto J3_btagCSVV2 = ROOT::VecOps::Take(Jet_btagCSVV2_masked, Trijet[2]);
            return J1_btagCSVV2 > 0.5 || J2_btagCSVV2 > 0.5 || J3_btagCSVV2 > 0.5;
            """
            # FIXME 
            # Do insteam something like max(J1_btag,J2_btag, Jt3_btag)>0.5. 
            # Do I need to define custom function max(RVec 1, RVec2, RVec3)?
        )

    # Assign four-momentums to each trijet combination
    df = (
        df.Define('J1', 'ROOT::VecOps::Take(Jet_p4, Trijet[0])').
        Define('J2', 'ROOT::VecOps::Take(Jet_p4, Trijet[1])').
        Define('J3', 'ROOT::VecOps::Take(Jet_p4, Trijet[2])').
        Define('Trijet_p4', '(J1+J2+J3)[Trijet_btag]')
    )



    # Get trijet transverse momentum values from four-momentum vectors
    df = df.Define(
        "Trijet_pt",
        "return ROOT::VecOps::Map(Trijet_p4, [](const ROOT::Math::PxPyPzMVector &v) { return v.Pt(); })",
    )

    # Evaluate mass of trijet with maximum pt and btag higher than threshold

    df = df.Define(
        "Trijet_mass", "Trijet_p4[ROOT::VecOps::ArgMax(Trijet_pt)].M()"
    )

    return df


def book_histos(
    df: ROOT.RDataFrame,
    process: str,
    variation: str,
    nevents: int,
) -> list[AGCResult]:
    """Return the RDataFrame results pertaining to the desired process and variation."""
    # Calculate normalization for MC
    x_sec = XSEC_INFO[process]
    lumi = 3378  # /pb
    xsec_weight = x_sec * lumi / nevents
    df = df.Define("Weights", str(xsec_weight))  # default weights

    if variation == "nominal":
        # Jet_pt variations definition
        # pt_scale_up() and pt_res_up(jet_pt) return scaling factors applying to jet_pt
        # pt_scale_up() - jet energy scaly systematic
        # pt_res_up(jet_pt) - jet resolution systematic
        df = df.Vary(
            "Jet_pt",
            "ROOT::RVec<ROOT::RVecF>{Jet_pt*pt_scale_up(), Jet_pt*jet_pt_resolution(Jet_pt.size())}",
            ["pt_scale_up", "pt_res_up"],
        )

        if process == "wjets":
            # Flat weight variation definition
            df = df.Vary(
                "Weights",
                "Weights*flat_variation()",
                [f"scale_var_{direction}" for direction in ["up", "down"]],
            )

    # Event selection - the core part of the algorithm applied for both regions
    # Selecting events containing at least one lepton and four jets with pT > 25 GeV
    # Applying requirement at least one of them must be b-tagged jet (see details in the specification)
    df = (
        df.Define("Electron_pt_mask", "Electron_pt>25")
        .Define("Muon_pt_mask", "Muon_pt>25")
        .Define("Jet_pt_mask", "Jet_pt>25")
        .Filter("Sum(Electron_pt_mask) + Sum(Muon_pt_mask) == 1")
        .Filter("Sum(Jet_pt_mask) >= 4")
    )

    # b-tagging variations for nominal samples
    if variation == "nominal":
        df = df.Vary(
            "Weights",
            # FIXME: Jet_pt[Jet_pt_mask] is called 4 times
            "ROOT::RVecD{Weights*btag_weight_variation(Jet_pt[Jet_pt_mask])}",
            [
                f"{weight_name}_{direction}"
                for weight_name in [f"btag_var_{i}" for i in range(4)]
                for direction in ["up", "down"]
            ],
        )

    # Define HT observable for the 4j1b region
    # Only one b-tagged region required
    # The observable is the total transvesre momentum
    # fmt: off

    # not strict condition is used because the same selection cut is applied in the reference implementation
    # https://github.com/iris-hep/analysis-grand-challenge/blob/main/analyses/cms-open-data-ttbar/ttbar_analysis_pipeline.py#L254
    # FIXME: Jet_btagCSVV2[Jet_pt_mask] is called 3 times
    df4j1b = df.Filter("Sum(Jet_btagCSVV2[Jet_pt_mask]>=0.5)==1")\
               .Define("HT", "Sum(Jet_pt[Jet_pt_mask])")
    # fmt: on

    # Define trijet_mass observable for the 4j2b region (this one is more complicated)
    df4j2b = define_trijet_mass(df)

    # Select the right VariationsFor function depending on RDF or DistRDF
    if type(df).__module__ == "DistRDF.Proxy":
        variationsfor_func = ROOT.RDF.Experimental.Distributed.VariationsFor
    else:
        variationsfor_func = ROOT.RDF.Experimental.VariationsFor

    # Book histograms and, if needed, their systematic variations
    results = []
    for df, observable, region in zip([df4j1b, df4j2b], ["HT", "Trijet_mass"], ["4j1b", "4j2b"]):
        histo_model = ROOT.RDF.TH1DModel(
            name=f"{region}_{process}_{variation}", title=process, nbinsx=25, xlow=50, xup=550
        )
        nominal_histo = df.Histo1D(histo_model, observable, "Weights")

        if variation == "nominal":
            varied_histos = variationsfor_func(nominal_histo)
            results.append(AGCResult(varied_histos, region, process, variation, nominal_histo))
        else:
            results.append(AGCResult(nominal_histo, region, process, variation, nominal_histo))
        print(f"Booked histogram {histo_model.fName}")

    # Return the booked results
    # Note that no event loop has run yet at this point (RDataFrame is lazy)
    return results


def load_cpp():
    """Load C++ helper functions. Works for both local and distributed execution."""
    try:
        # when using distributed RDataFrame 'helpers.cpp' is copied to the local_directory
        # of every worker (via `distribute_unique_paths`)
        localdir = get_worker().local_directory
        cpp_source = Path(localdir) / "helpers.cpp"
    except ValueError:
        # must be local execution
        cpp_source = "helpers.cpp"

    ROOT.gSystem.CompileMacro(str(cpp_source), "kO")


def main() -> None:
    program_start = time()
    args = parse_args()

    # Do not add histograms to TDirectories automatically: we'll do it ourselves as needed.
    ROOT.TH1.AddDirectory(False)
    # Disable interactive graphics: avoids canvases flashing on screen before we save them to file
    ROOT.gROOT.SetBatch(True)

    if args.verbose:
        # Set higher RDF verbosity for the rest of the program.
        # To only change the verbosity in a given scope, use ROOT.Experimental.RLogScopedVerbosity.
        ROOT.Detail.RDF.RDFLogChannel.SetVerbosity(ROOT.Experimental.ELogLevel.kInfo)

    if args.scheduler == "mt":
        # Setup for local, multi-thread RDataFrame
        ROOT.EnableImplicitMT(args.ncores)
        print(f"Number of threads: {ROOT.GetThreadPoolSize()}")
        client = None
        load_cpp()
        run_graphs = ROOT.RDF.RunGraphs
    else:
        # Setup for distributed RDataFrame
        client = create_dask_client(args.scheduler, args.ncores, args.hosts)
        ROOT.RDF.Experimental.Distributed.initialize(load_cpp)
        run_graphs = ROOT.RDF.Experimental.Distributed.RunGraphs

    # Book RDataFrame results
    inputs: list[AGCInput] = retrieve_inputs(
        args.n_max_files_per_sample, args.remote_data_prefix, args.data_cache
    )
    results: list[AGCResult] = []
    for input in inputs:
        df = make_rdf(input.paths, client, args.npartitions)
        results += book_histos(df, input.process, input.variation, input.nevents)
    print(f"Building the computation graphs took {time() - program_start:.2f} seconds")

    # Run the event loops for all processes and variations here
    run_graphs_start = time()
    run_graphs([r.nominal_histo for r in results])
    print(f"Executing the computation graphs took {time() - run_graphs_start:.2f} seconds")
    if client is not None:
        client.close()

    results = postprocess_results(results)
    save_plots(results)
    save_histos([r.histo for r in results], output_fname=args.output)
    print(f"Result histograms saved in file {args.output}")


if __name__ == "__main__":
    main()

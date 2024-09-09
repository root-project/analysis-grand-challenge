import argparse
import multiprocessing
from pathlib import Path
from time import time
from typing import Tuple

import ml
import ROOT
from distributed import Client, LocalCluster, SSHCluster, get_worker
from plotting import save_ml_plots, save_plots
from statistical import fit_histograms
from utils import AGCInput, AGCResult, postprocess_results, retrieve_inputs, save_histos

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
                `--remote-data-prefix='root://eospublic.cern.ch//eos/root-eos/AGC'`.""",
    )
    p.add_argument(
        "--output",
        "-o",
        help="Name of the file where analysis results will be stored. If it already exists, contents are overwritten.",
        default="histograms.root",
    )
    p.add_argument(
        "--inference",
        action=argparse.BooleanOptionalAction,
        help="""Produce machine learning histograms if enabled.
                        Disabled by default.""",
    )
    p.add_argument(
        "--scheduler",
        "-s",
        help="""The scheduler for RDataFrame parallelization. Will honor the --ncores flag.
                The default is `mt`, i.e. multi-thread execution.
                If dask-ssh, a list of worker node hostnames to connect to should be provided via the --nodes option.""",
        default="mt",
        choices=["mt", "dask-local", "dask-ssh", "dask-remote"],
    )
    p.add_argument(
        "--scheduler-address",
        help="Full address of the Dask scheduler, passed as argument to the distributed.Client object. If this argument is provided, the 'dask-remote' scheduler must be chosen, and it is a required argument in such case.",
    )
    p.add_argument(
        "--ncores",
        "-c",
        help=(
            "Number of cores to use. In case of distributed execution this is the amount of cores per node."
        ),
        default=multiprocessing.cpu_count(),
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
    p.add_argument(
        "-v",
        "--verbose",
        help="Turn on verbose execution logs.",
        action="store_true",
    )

    p.add_argument(
        "--statistical-validation",
        help = argparse.SUPPRESS,
        action="store_true",
    )

    p.add_argument(
        "--no-fitting",
        help="Do not run statistical validation part of the analysis.",
        action="store_true",
    )

    return p.parse_args()


def create_dask_client(scheduler: str, ncores: int, hosts: str, scheduler_address: str) -> Client:
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
            worker_options={
                "nprocs": ncores,
                "nthreads": 1,
                "memory_limit": "32GB",
            },
        )
        return Client(sshc)

    if scheduler == "dask-remote":
        return Client(scheduler_address)

    raise ValueError(
        f"Unexpected scheduling mode '{scheduler}'. Valid modes are ['dask-local', 'dask-ssh', 'dask-remote']."
    )


def define_trijet_mass(df: ROOT.RDataFrame) -> ROOT.RDataFrame:
    """Add the trijet_mass observable to the dataframe after applying the appropriate selections."""

    # First, select events with at least 2 b-tagged jets
    df = df.Filter("Sum(Jet_btagCSVV2_cut > 0.5) > 1")

    # Build four-momentum vectors for each jet
    df = df.Define(
        "Jet_p4",
        "ConstructP4(Jet_pt_cut, Jet_eta_cut, Jet_phi_cut, Jet_mass_cut)",
    )

    # Build trijet combinations
    df = df.Define("Trijet_idx", "Combinations(Jet_pt_cut, 3)")

    # Trijet_btag is a helpful array mask indicating whether or not the maximum btag value in Trijet is larger than the 0.5 threshold
    df = df.Define(
        "Trijet_btag",
        """
            auto J1_btagCSVV2 = Take(Jet_btagCSVV2_cut, Trijet_idx[0]);
            auto J2_btagCSVV2 = Take(Jet_btagCSVV2_cut, Trijet_idx[1]);
            auto J3_btagCSVV2 = Take(Jet_btagCSVV2_cut, Trijet_idx[2]);
            return J1_btagCSVV2 > 0.5 || J2_btagCSVV2 > 0.5 || J3_btagCSVV2 > 0.5;
            """,
    )

    # Assign four-momentums to each trijet combination
    df = df.Define(
        "Trijet_p4",
        """
        auto J1 = Take(Jet_p4, Trijet_idx[0]);
        auto J2 = Take(Jet_p4, Trijet_idx[1]);
        auto J3 = Take(Jet_p4, Trijet_idx[2]);
        return (J1+J2+J3)[Trijet_btag];
        """,
    )

    # Get trijet transverse momentum values from four-momentum vectors
    df = df.Define(
        "Trijet_pt",
        "return Map(Trijet_p4, [](const ROOT::Math::PxPyPzMVector &v) { return v.Pt(); })",
    )

    # Evaluate mass of trijet with maximum pt and btag higher than threshold
    df = df.Define("Trijet_mass", "Trijet_p4[ArgMax(Trijet_pt)].M()")

    return df


def book_histos(
    df: ROOT.RDataFrame, process: str, variation: str, nevents: int, inference=False
) -> Tuple[list[AGCResult], list[AGCResult]]:
    """Return the pair of lists of RDataFrame results pertaining to the given process and variation.
    The first list contains histograms of reconstructed HT and trijet masses.
    The second contains ML inference outputs"""
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
            "ROOT::RVec<ROOT::RVecF>{Jet_pt*pt_scale_up(), Jet_pt*jet_pt_resolution(Jet_pt)}",
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
    # Selecting events containing at least one lepton and four jets with pT > 30 GeV
    # Applying requirement at least one of them must be b-tagged jet (see details in the specification)
    df = (
        df.Define(
            "Electron_mask",
            "Electron_pt > 30 && abs(Electron_eta) < 2.1 && Electron_sip3d < 4 && Electron_cutBased == 4",
        )
        .Define(
            "Muon_mask",
            "Muon_pt > 30 && abs(Muon_eta) < 2.1 && Muon_sip3d < 4 && Muon_tightId && Muon_pfRelIso04_all < 0.15",
        )
        .Filter("Sum(Electron_mask) + Sum(Muon_mask) == 1")
        .Define("Jet_mask", "Jet_pt > 30 && abs(Jet_eta) < 2.4 && Jet_jetId == 6")
        .Filter("Sum(Jet_mask) >= 4")
    )

    # create columns for "good" jets
    df = (
        df.Define("Jet_pt_cut", "Jet_pt[Jet_mask]")
        .Define("Jet_btagCSVV2_cut", "Jet_btagCSVV2[Jet_mask]")
        .Define("Jet_eta_cut", "Jet_eta[Jet_mask]")
        .Define("Jet_phi_cut", "Jet_phi[Jet_mask]")
        .Define("Jet_mass_cut", "Jet_mass[Jet_mask]")
    )

    # b-tagging variations for nominal samples
    if variation == "nominal":
        df = df.Vary(
            "Weights",
            "ROOT::RVecD{Weights*btag_weight_variation(Jet_pt_cut)}",
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
    df4j1b = df.Filter("Sum(Jet_btagCSVV2_cut > 0.5) == 1").Define("HT", "Sum(Jet_pt_cut)")
    # fmt: on

    # Define trijet_mass observable for the 4j2b region (this one is more complicated)
    df4j2b = define_trijet_mass(df)

    # Book histograms and, if needed, their systematic variations
    results = []
    for df, observable, region in zip([df4j1b, df4j2b], ["HT", "Trijet_mass"], ["4j1b", "4j2b"]):
        histo_model = ROOT.RDF.TH1DModel(
            name=f"{region}_{process}_{variation}",
            title=process,
            nbinsx=25,
            xlow=50,
            xup=550,
        )
        nominal_histo = df.Histo1D(histo_model, observable, "Weights")

        if variation == "nominal":
            results.append(
                AGCResult(
                    nominal_histo,
                    region,
                    process,
                    variation,
                    nominal_histo,
                    should_vary=True,
                )
            )
        else:
            results.append(
                AGCResult(
                    nominal_histo,
                    region,
                    process,
                    variation,
                    nominal_histo,
                    should_vary=False,
                )
            )
        print(f"Booked histogram {histo_model.fName}")

    ml_results: list[AGCResult] = []

    if not inference:
        return (results, ml_results)

    df4j2b = ml.define_features(df4j2b)
    df4j2b = ml.infer_output_ml_features(df4j2b)

    # Book histograms and, if needed, their systematic variations
    for i, feature in enumerate(ml.ml_features_config):
        histo_model = ROOT.RDF.TH1DModel(
            name=f"{feature.name}_{process}_{variation}",
            title=feature.title,
            nbinsx=feature.binning[0],
            xlow=feature.binning[1],
            xup=feature.binning[2],
        )

        nominal_histo = df4j2b.Histo1D(histo_model, f"results{i}", "Weights")

        if variation == "nominal":
            ml_results.append(
                AGCResult(
                    nominal_histo,
                    feature.name,
                    process,
                    variation,
                    nominal_histo,
                    should_vary=True,
                )
            )
        else:
            ml_results.append(
                AGCResult(
                    nominal_histo,
                    feature.name,
                    process,
                    variation,
                    nominal_histo,
                    should_vary=False,
                )
            )
        print(f"Booked histogram {histo_model.fName}")

    # Return the booked results
    # Note that no event loop has run yet at this point (RDataFrame is lazy)
    return (results, ml_results)


def load_cpp():
    """Load C++ helper functions. Works for both local and distributed execution."""
    try:
        # when using distributed RDataFrame 'helpers.cpp' is copied to the local_directory
        # of every worker (via `distribute_unique_paths`)
        localdir = get_worker().local_directory
        cpp_source = Path(localdir) / "helpers.h"
    except ValueError:
        # must be local execution
        cpp_source = "helpers.h"

    ROOT.gInterpreter.Declare(f'#include "{str(cpp_source)}"')


def run_mt(
    program_start: float,
    args: argparse.Namespace,
    inputs: list[AGCInput],
    results: list[AGCResult],
    ml_results: list[AGCResult],
) -> None:
    ROOT.EnableImplicitMT(args.ncores)
    print(f"Number of threads: {ROOT.GetThreadPoolSize()}")
    load_cpp()
    if args.inference:
        ml.load_cpp()

    for input in inputs:
        df = ROOT.RDataFrame("Events", input.paths)
        hist_list, ml_hist_list = book_histos(
            df, input.process, input.variation, input.nevents, inference=args.inference
        )
        results += hist_list
        ml_results += ml_hist_list

    for r in results + ml_results:
        if r.should_vary:
            r.histo = ROOT.RDF.Experimental.VariationsFor(r.histo)

    print(f"Building the computation graphs took {time() - program_start:.2f} seconds")

    # Run the event loops for all processes and variations here
    run_graphs_start = time()
    ROOT.RDF.RunGraphs([r.nominal_histo for r in results + ml_results])

    print(f"Executing the computation graphs took {time() - run_graphs_start:.2f} seconds")


def run_distributed(
    program_start: float,
    args: argparse.Namespace,
    inputs: list[AGCInput],
    results: list[AGCResult],
    ml_results: list[AGCResult],
) -> None:
    if args.inference:

        def ml_init():
            load_cpp()
            ml.load_cpp()

        ROOT.RDF.Experimental.Distributed.initialize(ml_init)
    else:
        ROOT.RDF.Experimental.Distributed.initialize(load_cpp)

    scheduler_address = args.scheduler_address if args.scheduler_address else ""
    with create_dask_client(args.scheduler, args.ncores, args.hosts, scheduler_address) as client:
        for input in inputs:
            df = ROOT.RDF.Experimental.Distributed.Dask.RDataFrame(
                "Events",
                input.paths,
                daskclient=client,
                npartitions=args.npartitions,
            )
            df._headnode.backend.distribute_unique_paths(
                [
                    "helpers.h",
                    "ml_helpers.h",
                    "ml.py",
                    "models/bdt_even.root",
                    "models/bdt_odd.root",
                ]
            )
            hist_list, ml_hist_list = book_histos(
                df, input.process, input.variation, input.nevents, inference=args.inference
            )
            results += hist_list
            ml_results += ml_hist_list

        for r in results + ml_results:
            if r.should_vary:
                r.histo = ROOT.RDF.Experimental.Distributed.VariationsFor(r.histo)

        print(f"Building the computation graphs took {time() - program_start:.2f} seconds")

        # Run the event loops for all processes and variations here
        run_graphs_start = time()
        ROOT.RDF.Experimental.Distributed.RunGraphs([r.nominal_histo for r in results + ml_results])

        print(f"Executing the computation graphs took {time() - run_graphs_start:.2f} seconds")


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
        ROOT.Detail.RDF.RDFLogChannel().SetVerbosity(ROOT.Experimental.ELogLevel.kInfo)

    if args.statistical_validation:
        fit_histograms(filename=args.output)
        return

    inputs: list[AGCInput] = retrieve_inputs(
        args.n_max_files_per_sample, args.remote_data_prefix, args.data_cache
    )
    results: list[AGCResult] = []
    ml_results: list[AGCResult] = []

    if args.scheduler == "mt":
        run_mt(program_start, args, inputs, results, ml_results)
    else:
        if args.scheduler == "dask-remote" and not args.scheduler_address:
            raise ValueError(
                "'dask-remote' option chosen but no address provided for the scheduler. Provide it with `--scheduler-address`."
            )
        if args.scheduler_address and args.scheduler != "dask-remote":
            raise ValueError(
                f"An address of a Dask scheduler was provided but the chosen scheduler is '{args.scheduler}'. The 'dask-remote' scheduler must be chosen if an address is provided."
            )
        run_distributed(program_start, args, inputs, results, ml_results)

    results = postprocess_results(results)
    save_plots(results)
    save_histos([r.histo for r in results], output_fname=args.output)
    print(f"Result histograms saved in file {args.output}")

    if args.inference:
        ml_results = postprocess_results(ml_results)
        save_ml_plots(ml_results)
        output_fname = args.output.split(".root")[0] + "_ml_inference.root"
        save_histos([r.histo for r in ml_results], output_fname=output_fname)
        print(f"Result histograms from ML inference step saved in file {output_fname}")

    if not args.no_fitting:
        fit_histograms(filename=args.output)


if __name__ == "__main__":
    main()

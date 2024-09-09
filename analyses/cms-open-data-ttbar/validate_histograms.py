# This script is originally from the IRIS-HEP analysis-grand-challenge repository
# and has been adapted to serve the purpose of the RDF-based analysis.
# The original script can be found here:
# https://github.com/iris-hep/analysis-grand-challenge/blob/main/analyses/cms-open-data-ttbar/validate_histograms.py

# Compare the content of histograms produced by analysis with a reference file.
# A reference file for N_FILES_MAX_PER_SAMPLE=1 is available in directory `reference/`.

from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict

import numpy as np
import uproot


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--histos",
        help="ROOT file containing the output histograms. Defaults to './histograms.root'.",
        default="histograms.root",
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--reference",
        help="JSON reference against which histogram contents should be compared",
    )
    group.add_argument(
        "--dump-json",
        help="Print JSON representation of histogram contents to screen",
        action="store_true",
    )
    parser.add_argument(
        "--verbose",
        help="Print extra information about bin mismatches",
        action="store_true",
    )
    return parser.parse_args()


# convert uproot file containing only TH1Ds to a corresponding JSON-compatible dict with structure:
# { "histo1": { "edges": [...], "contents": [...] }, "histo2": { ... }, ... }
# Only the highest namecycle for every histogram is considered, and cycles are stripped from the histogram names.
def as_dict(f: uproot.ReadOnlyDirectory) -> dict[str, dict]:
    histos: dict[str, dict] = defaultdict(dict)
    # this assumes that the rightmost ";" (if any) comes before a namecycle
    names = set(k.rsplit(";", 1)[0] for k in f)
    for name in names:
        h = f[name]
        assert isinstance(h, uproot.behaviors.TH1.Histogram)
        histos[name]["edges"] = h.axis().edges().tolist()
        histos[name]["contents"] = h.counts(flow=True).tolist()
    return histos


def find_matching_histogram(name: str, histos: dict) -> str | None:
    if name in histos:
        return name
    for hname in histos.keys():
        if hname.startswith(name):
            return hname
    return None


def validate(
    histos: dict, reference: dict, verbose=False
) -> dict[str, list[str]]:
    errors = defaultdict(list)
    for name, ref_h in reference.items():
        matching_name = find_matching_histogram(name, histos)
        if matching_name is None:
            if name != "4j1b_pseudodata" and name != "4j2b_pseudodata":
                errors[name].append("Histogram not found.")
            continue

        h = histos[matching_name]
        if not np.allclose(h["edges"], ref_h["edges"]):
            errors[name].append(
                f"Edges do not match:\n\tgot      {h['edges']}\n\texpected {ref_h['edges']}"
            )
        contents_depend_on_rng = "pt_res_up" in name
        is_close = np.isclose(h["contents"], ref_h["contents"])

        if not contents_depend_on_rng and not all(is_close):
            where_not_close = np.where(np.invert(is_close))[0]
            diff_values = np.diff(where_not_close)
            split_indices = np.argwhere(np.abs(diff_values) > 1).reshape(-1)

            if len(split_indices) == 0:
                split_values = [where_not_close]
            else:
                split_indices += 1
                split_values = np.split(where_not_close, split_indices)

            is_error = False
            for group in split_values:
                h_group = np.array(h["contents"])[group]
                ref_group = np.array(ref_h["contents"])[group]
                if not np.allclose(h_group, ref_group, atol=2.0):
                    is_error = True
                    if verbose:
                        print(f"In {name}: Not close enough for bin migration")
                        print(
                            "histogram: ", h_group, ", reference: ", ref_group
                        )
                        print()
                elif not np.allclose(sum(h_group), sum(ref_group)):
                    is_error = True
                    if verbose:
                        print(f"In {name}: Partial sums are unequal")
                        print(
                            "histogram: ", h_group, ", reference: ", ref_group
                        )
                        print()
                else:
                    if verbose:
                        print(f"In {name}: Bin migration likely")
                        print(
                            "histogram: ", h_group, ", reference: ", ref_group
                        )
                        print()
            if is_error:
                errors[name].append(
                    f"Contents do not match:\n\tgot      {h['contents']}\n\texpected {ref_h['contents']}"
                )

    return errors


if __name__ == "__main__":
    args = parse_args()
    with uproot.open(args.histos) as f:
        histos = as_dict(f)

    if args.dump_json:
        print(json.dumps(histos, indent=2, sort_keys=True))
        sys.exit(0)

    with open(args.reference) as reference:
        ref_histos = json.load(reference)

    print(
        f"Validating '{args.histos}' against reference '{args.reference}'..."
    )
    errs = validate(histos=histos, reference=ref_histos, verbose=args.verbose)

    if len(errs) == 0:
        print("All good!")
    else:
        for hist_name, errors in errs.items():
            errstr = '\n\t'.join(errors)
            print(f"{hist_name}\n\t{errstr}")
        sys.exit(1)

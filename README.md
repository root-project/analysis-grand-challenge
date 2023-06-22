# Analysis Grand Challenge (AGC) benchmarks with ROOT

The Analysis Grand Challenge (AGC) is about performing the last steps in an analysis pipeline at scale to test workflows envisioned for the HL-LHC.
This includes:

- columnar data extraction from large datasets
- processing of that data (event filtering, construction of observables, evaluation of systematic uncertainties) into histograms
- statistical model construction and statistical inference
- visualizations for these steps

all done in a reproducible & preservable way that can scale to HL-LHC requirements.

The official AGC documentation can be found [at this link](https://agc.readthedocs.io/en/latest/index.html).

This repository implements AGC analysis tasks using modern ROOT interfaces.

The physics analyses implementations can be found in the `analyses` folder.

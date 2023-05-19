# Analysis Grand Challenge (AGC) benchmarks with ROOT

The Analysis Grand Challenge (AGC) is about performing the last steps in an analysis pipeline at scale to test workflows envisioned for the HL-LHC.
This includes

- columnar data extraction from large datasets,
- processing of that data (event filtering, construction of observables, evaluation of systematic uncertainties) into histograms,
- statistical model construction and statistical inference,
- relevant visualizations for these steps,

all done in a reproducible & preservable way that can scale to HL-LHC requirements.

The reference implementation of the workflows is hosted at https://github.com/iris-hep/analysis-grand-challenge. This repository demonstrates usage of ROOT facilities in the same type of benchmarks. The physics analyses code can be found in the `analyses` folder.

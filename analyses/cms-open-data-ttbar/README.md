The physics analysis task is a $t\bar{t}$ cross-section measurement with 2015 CMS Open Data. Currently, the benchmark is fixed at the reference implementation tag [v0.1.0](https://github.com/iris-hep/analysis-grand-challenge/tree/v0.1.0). The first RDataFrame implementation of this task is a work by Andrii Falko and can be found at https://github.com/andriiknu/analysis-grand-challenge.

The benchmark assumes the latest ROOT version 6.28 is used. For installation instructions, refer to https://root.cern/install.

The sub-folders represent different modes of execution:

- `rdf-distributed`: run the analysis on one or multiple computing nodes with multiprocessing, leveraging the [Dask](https://www.dask.org/) distributed backend of RDataFrame.
- `rdf-imt`: run the analysis using the full power of a single machine through implicit multithreading parallelisation, leveraging the native C++ multithreading capabilities within ROOT.

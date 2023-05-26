The physics analysis task is a $t\bar{t}$ cross-section measurement with 2015 CMS Open Data. Currently, the benchmark is fixed at the reference implementation tag [v0.1.0](https://github.com/iris-hep/analysis-grand-challenge/tree/v0.1.0). The first RDataFrame implementation of this task is a work by Andrii Falko and can be found at https://github.com/andriiknu/analysis-grand-challenge.

The benchmark assumes the latest ROOT version 6.28 is used. For installation instructions, refer to https://root.cern/install.

To run the analysis:

```
python analysis.py [OPTIONS ...]
```

To see the full list of options, run `python analysis.py -h`. RDataFrame supports local or distributed mode, which can be toggled via `--scheduling-mode [imt, dask-local, ...]`.
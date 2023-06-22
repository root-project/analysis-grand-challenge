# $t\bar{t}$ analysis

The physics analysis task is a $t\bar{t}$ cross-section measurement with 2015 CMS Open Data.

A full explanation of the analysis task is available [at this link](https://agc.readthedocs.io/en/latest/taskbackground.html).

The current implementation corresponds to [AGC version 0.1.0](https://github.com/iris-hep/analysis-grand-challenge/tree/v0.1.0).

## Running the analysis

ROOT v6.28 or later is required (see https://root.cern/install for installation instructions).

The full analysis can be run in multi-thread mode with:

```
python analysis.py
```

Use `python analysis.py -h` to see the full list of options, including how to switch between local multi-threading and distributed execution.

## Acknowledgements

The first RDataFrame implementation of this task was a work by Andrii Falko (@andriiknu) sponsored by the IRIS-HEP Fellows Program.

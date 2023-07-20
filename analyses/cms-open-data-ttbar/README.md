# $t\bar{t}$ analysis

The physics analysis task is a $t\bar{t}$ cross-section measurement with 2015 CMS Open Data.

A full explanation of the analysis task is available [at this link](https://agc.readthedocs.io/en/latest/taskbackground.html).

## Versions

- [tag v1.0.0](https://github.com/root-project/analysis-grand-challenge/tree/v1.0.0) of this repository implements v1 of the AGC task and requires ROOT v6.28.04 or above; the histograms produced are expected to be bin-by-bin compatible with [v1.2.0](https://github.com/iris-hep/analysis-grand-challenge/tree/v1.2.0) of the reference implementation
- the [main branch](https://github.com/root-project/analysis-grand-challenge/tree/main) of this repository is a work in progress towards v2 of the AGC task; it requires a [ROOT nightly build](root.cern/nightlies)

More about [AGC versions](https://agc.readthedocs.io/en/latest/versionsdescription.html).

## Running the analysis

Given a working installation of the appropriate ROOT version (see above) the full analysis can be run in multi-thread mode with:

```
python analysis.py
```

Use `python analysis.py -h` to see the full list of options, including how to switch between local multi-threading and distributed execution and how to limit the execution to a subset of the full dataset.

## Acknowledgements

The first RDataFrame implementation of this task was a work by Andrii Falko (@andriiknu) sponsored by the IRIS-HEP Fellows Program.

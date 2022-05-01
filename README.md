# zebrafinchRA
Scripts to analyze developmental and sex-biased gene expression from zebra finch RA

## Getting started

The easiest way to get going is to load dependencies via the renv virtual environment included with this project. First install the renv package to your machine by running `packages.install("renv")` in RStudio. Then, from the main project folder (`zebrafinchRA/`), run `renv::restore()` to setup a virtual environment with all of the packages as declared in the lockfile (`renv.lock`). Finally, run the `deseq_analysis.R` file in steps or on the whole.

Note: This project relies on the xlsx package to write results tables to Excel files. The xlsx package has Java dependencies, so Java must be installed on the user's machine for these results table export functions to work. If this is not feasible, the R results table objects contain all data reported in [paper tables] and can be exported to other formats through alternative methods/packages (e.g. readxls).

## A note on reproducibility

The `deseq_analysis.R` script contains several plotting functions which produce the figures as presented in [paper].

The output tables for GO analyses may vary slightly from run to run due to random sampling, as well as changes to the GO mappings for examined genes. However, the overall themes reported in Figure XX of [paper] persisted across multiple runs.

## Raw data availability

Raw sequencing data can be found on the GEO ...
# Utilities for Y chromosome analyses

## Overview

This repository contains a collection of utility scripts and tools developed primarily for the paper: *Resolving the source of branch length variation in the Y chromosome phylogeny*.  The utilities involve Python scripts (`VcfHandler.py`, `BranchShortening.py`, `YChrDataset.py`, `BamUtils.py`) and an R script (`utils.R`). As these scripts are used elsewhere, this repo is designed to be a submodule that can be included along with the project-specific analyses.

`VcfHandler.py` contains a number of functions used to read and process VCF files (typically generated with snpAD) and .bed files.

`BranchShortening.py` contains functions used to calculate Y chromosome TMRCAs and analyse Y chromosome reference bias.

`YChrDataset.py` stores metadata about the Y chromosome dataset analysed in the paper.

`BamUtils.py` contains functions used to process .bam files.

`utils.R` stores information about the Y chromosome dataset and is used for plotting trees in R using `ggplot`.

## Installation
The required python packages are `pandas`, `numpy`, `statsmodels`, `matplotlib` and `pydantic`

Some functions in `BamUtils.py` require that `samtools` is installed and is included in `PATH`.
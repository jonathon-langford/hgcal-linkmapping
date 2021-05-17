# HGCal link mapping studies

## Introduction

This README details how to obtain the configuration of trigger links between the CMS HGCal and Stage 1 FPGAs.
A recipe to run each of the stages is given as well as some information on some useful auxiliary functions.
In each case the default config files show which functions and options are available.

## Installation

- If you don't have python and the relevant packages installed please run: `source install_packages.sh`.
Note the size of the installation is around 4GB.
Then each time you start a new session run: `source start_mapping_env.sh`

- If you have python and can install your own packages the relevant ones are:
    - numpy
    - matplotlib	
    - pandas
    - pyyaml
    - scikit-learn
    - ROOT
    - root_numpy

## Obtaining the input histograms

One of the inputs to the minimisation process (described later) is a `ROOT` file containing 2D histograms of `phi` vs `r/z` for each module in "sector 0" (i.e. one 120 degree slice) of one end-cap. The code doing this is in `extract_data.cxx`.

This is run in the following manner:
-  `make`
-  `./extract_data config/extract_data.json`

As input (`inputfile`) the program takes an output `ROOT` ntuple from the CMSSW TPG framework. The name of the output histogram file is given by `file_ROverZHistograms`. The binning of the 2D histograms in r/z is defined by the variables `nROverZBins`, `rOverZMin`, and `rOverZMax`.

There is the option to create 2D histograms suitable for use with the `FeMappingV7` or `FeMappingTpgV7` mapping file (`configFileVersion` parameter).

The other variables set the name of some auxiliary outputs (which are not needed for the main recipe).

## `main.py`

Main file containing the option to run all functions:
- `plot_lpGBTLoads`, `plot_ModuleLoads`, : Processes MC event data and determines the average number of TCs, or words, per lpGBT
- `study_mapping`, :  Find the optimised way of assigning lpGBTs to bundles

Run using the config file `config/default.yaml`, where the input options are listed for each parameter

## `process.py`

Contains the helper functions required for each function in `main.py`

## `rotate.py` and `rotate.cxx`

Python and C++ implementations of the mapping between 120 degree HGCal sectors in (u,v) coordinates.

## `fluctuation.py`

Takes as input a choice of lpgbt bundles, and bins the trigger cell data event by event
There are several plotting scripts that investigate the impact of truncation on the number of trigger cells.
Run using the config file `config/fluctuation.yaml`.
Also the option to save the sum of (truncated or total) trigger cell p<sub>T</sub> as a function of R/Z for each event.

## `plotbundles.py`

Various plotting functions, mainly to plot the 24 R/Z histograms for each bundle, and take the ratio to the inclusive distribution over 24.
Run using the config file `config/plotbundles.yaml`
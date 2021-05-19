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

## Performing the minimisation

The function relevant for performing the minimisation are found in `main.py` (which imports and uses additional functions from `process.py`.

The code is run in the following manner:
- ./main.py config/default.yaml i

where `i` is a number which is attached to the output configuration on completion. This is useful when running many parallel instances on a batch system.

There are seven possible functions in `main.py` and these are run by setting the relevant function to `True` in `config/default.py`. The most important function is `study_mapping`. This needs to be `True` in order to perform the minimisation and find the optimised way of assigning lpGBTs to FPGAs (also known as bundles). The other functions are described at the end of this README.

The additional configurable parameters for the `study_mapping` function are in the `study_mapping` block of the config file.

`MappingFile` gives the input location for the file listing each module and which lpGBTs are connected to each. `CMSSW_ModuleHists` gives the input location for the set of 2D histograms that were created in the first step of `extract_data.cxx`. `TowerMappingFile` gives the location of the file listing each module and which towers overlap with each. This is used if `include_max_towers_in_chi2` is set to `True` below.

The minimisation is either performed using a random hill climb or simulated annealing algorithm. The choice is set under `algorithm`. Note that there is also an option `save_root`, which directly saves the r/z histograms for each bundle to a `ROOT` file (for the `initial_state` configuration). This option is not often used.

The initial state is either set to be `random` (in which case there is also the option to set the `random_seed` - otherwise `random_seed = ~` which means it is not set), or it is set to an initial configuration from an input file. The default file `data/bundles_example-pickle-tpgv7.npy` is provided. Any output file from the minimisation can be used as an input to another minimisation.

`max_iterations` defines how many iterations should be performed in the minimisation before ending. Note that a best-so-far configuration is saved to a file if the minimisation is ended before reaching a minimum (either by keyboard interrupt, or reaching the maximum number of iterations. The `minigroup_type` parameter defines the philosophy of forming the mini-groups, which are small groups of modules which must be treated together as one in the minimisation. Generally `minimal` should be used and is the most tested.

One has several options in defining the &Chi^2 used in the minimisation. One can use only the r/z values of each bundle histogram, rather than the associate statistical uncertainties (set `include_errors_in_chi2` to `False`).
By default the maximum number of modules attached to an FPGA is used in the &Chi^2 (`include_max_modules_in_chi2`) as well as the maximum number of towers covered by an FPGA (`include_max_towers_in_chi2`).
These values are included with an arbitrary constant weighting term (`max_modules_weighting_factor` and `max_towers_weighting_factor` respectively).

For the towers, there is an additional option `max_towers_weighting_option` which is set to 1 if the maximum is used directly in the &Chi^2 and 2 if it is used with a step-like penalty function. If `TowerPhiSplit` is not None or `[]` then the maximum number of towers covered by an FPGA in a certain phi region is the relevant quantity. The default is `[6,15]`. The number refers to 5 degree bins in phi (where 0-5 degrees is bin 1). [6,16] therefore means 3 phi regions, phi < 30 degrees, 30 < phi < 80 and phi > 80 degrees.

Finally there is an option `weight_bins_proportionally` which divides the &Chi^2 in each r/z bin by the bin value (multiplied by a constant).

The configurable information in `phisplit` details how to split the 2D r/z histograms in phi, so as to define the phidivisionX and phidivisionY regions. The default is `per_rover_bin`, which means the mid-point in phi in each r/z bin is used as the division. The other option is `fixed`, which means the split is at a fixed point in `phi` (the values of which are defined using the `phidivisionX_fixvalue_min` and `phidivisionX_fixvalue_max` variables.

Finally any corrections to account for differences between the geometry in the input `ROOT` histograms and the latest geometry are given in the `corrections` block. If using `v11` geometry these should generally be left unchanged.

## Plotting the best mapping

Once an output mapping configuration file has been obtained (with the `study_mapping` function described above), one might wish to plot the r/z histograms of the 25 bundles (FPGAs), and take the ratio to the inclusive distribution divided by 24. This is achieved using the `plotbundles.py` file and run like:

- `./plotbundles.py config/plotbundles.yaml`

The configuration file uses very similar options to the `default` configuration file. The main differences are the `input_file` (the output mapping configuration from `study_mapping`) and the `output_dir` (the location where the plots should be saved.

## Calculating truncation values and determining the effect of truncation, event by event.

### Filling the r/z histograms

Once a mapping configuration from the `study_mapping` function has been found, the next step is to estimate the amount of truncation that would be needed in each r/z bin in order to stay within the limits given by the hardware.
This is done using the two files `fluctuation.py` and `fluctuation_postprocess.py`. The first, `fluctuation.py` contains the `main` function as well as a function `checkFluctuations` to fill r/z histograms for each bundle in each event (where in practise the histograms are filled 6 times per event due to the six fold symmetry). 
This is run in the following manner:

- `./fluctuation.py config/fluctuations.yaml`

The relevant variables in the config file are now described.
As in the file `default.yaml` the chosen function is set to be `True` or `False` in the first block. In this case `checkFluctuations` should be set to `True`.
Some of the parameters are the same as in previous configs and are not described again.

The block `checkFluctuations` provides the additional configurable parameters for this function.
`beginEvent` and `endEvent` define the start and end events for running over.
The relevant input file is the same TPG CMSSW `ROOT` file that was also used to produce the r/z histograms from `extract_data.cxx`. In addition one requires the `initial_state`, which is the output mapping configuration from the minimisation and the `mappingFile` used in the minimisation. The `outputName` name is also set here. The parameter `save_ntc_hists` determines if auxiliary histograms showing the number of trigger cells per module are saved, by default this is `False`.

In addition to the familiar `phisplit` config there is also the `tcPtConfig`.
There is the option (if `save_sum_tcPt` = `True`) to save the sum of (truncated or total) trigger cell p<sub>T</sub> as a function of r/z for each event.
In this case one also needs to fill the  block `truncationConfig`, detailing the trigger cell limits and definition of the links.
Various options of interest can be defined. Each must give the maximum number of TCs on a link in regionA and regionB (`maxTCsA` and `maxTCsB` respectively). Then the number of links between Stage 1 and Stage 2 (`nLinks`) as well as how regionA and regionB are defined in terms of `phiDivisionX` and `phiDivisionY` (`regionADefinition` and `regionBDefinition`).
Finally it is necessary to give the TC truncation values in each r/z bin (as a list). These must be determined in an initial run (see below for how) as otherwise the output files are too large.

### Determining the truncation values

In order to determine the optimal trunction values one uses the function `studyTruncationOptions` in `fluctuation_postprocess.py`

It is run in the same manner as before, ensuring `studyTruncationOptions` is set to `True` in the config file:
- `./fluctuation.py config/fluctuations.yaml`

The relevant parameters are:
-  `eventData`, which is the output filled r/z histograms from the previous step
-  `options_to_study`, the numbers of the different scenarios under investigation. These numbers correspond to the option number in the `truncationConfig` described previously.
-  `truncation_values_method` - the two possibilities are `original` and `reverse`. The original method calculates the 99th percentile of TCs in each r/z bin and finds a scaling factor to apply to this to ensure that the sum of TCs is less than the maximum allowed. There is then a small redistribution of TCs. The reverse method aims to maintain a constant fraction of TCs truncated across r/z, and finds the highest possible fraction such that the sum of TCs is less than the maximum allowed. There is again a small redistribution of TCs.


## Further additional useful functions

In `main.py`:

- `plot_lpGBTLoads`, `plot_ModuleLoads`, : Processes MC event data and determines the average number of TCs, or words, per lpGBT

In `rotate.py` and `rotate.cxx`
- Python and C++ implementations of the mapping between 120 degree HGCal sectors in (u,v) coordinates.

In `fluctuation_postprocess.py`
- plot_MeanMax, 


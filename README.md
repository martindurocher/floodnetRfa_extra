# floodnetRfa extra

# Introduction

In this repository, I provide supplemental materials that add to my work on the R-package [floodnetRfa](https://github.com/floodnetProject16). 
It includes codes and summaries resulting from the analysis of the 1114 Canadian stations.
In particular, it includes part of the information I collected during the development of the R-package, which I summarize in two datasets: `gauged_sites` and `descriptors`.

The objective is to give further insights into the package functionalities and development. 
I derive all results from a version of the HYDAT database dating of August 11th, 2019. 
This specific version can be download [here](https://drive.google.com/file/d/1YI8pmB0U2Tp9FVVPpu2So8SmWIid9PsP/view?usp=sharing), and the user must add the local path of the database into the file `config` of the root folder.

The repository is divided into independent folders that I describe individually.
Note that the results are provided "as is",  and come without warranty of any kind. 

# Folders

## Autopot

One method available to perform at-site flood frequency analysis is Peaks Over Threshold (POT). 
The objective of this section is to use the automatic selection method provided in `floodnetRfa` to obtain candidates thresholds for each station in `gauged_sites.csv`. 
The script `autopot_script.R` selects the thresholds, and the file `autopot_readme.html` briefly summarizes them.

## Cache

In several sections, scripts run parallel tasks to speed up computation. 
These scripts save temporary results for each site in its cache folder. 
For instance, `RFA/RFA_script.R` uses the `cache/RFA` folder. 
The respective sections provide further details.

## Dataset

The folder contains two CSV files (zipped) that summarize the information collected on the 1114 Canadian hydrometric stations.

## Other

The section contains various small scripts. 

* The file `cmd_coordinates.R` investigates the use of Classical Multidimensional Scaling (CMD) to obtain a coordinate system that preserves the great-circle distance among studied sites. 

* The file `GetSuperRegions.R` contains a function that returns the super regions of a target site as suggested in the dataset `gauged_sites.csv`.

* The file `intersite.R` studies ways of estimating the coefficient of correlation of multivariate normal distribution that characterize the spatial dependence among annual maxima.

* The file `LmomDiag_fitting.R` find polynomial approximations of curves use in the L-moments ratio diagram. 

## RFA

At-site and regional frequency analysis are performed on each of the 1114 stations. 
The files `atsite_script.R` and `RFA_script.R` fit the extreme distributions and  `RFA_readme.Rmd` compares the precision of the estimators based on the confidence intervals lengths. 

## Super regions

The recommended strategy for forming pooling groups first identifies a set of stations having similar hydrological properties called super regions. The file `super_regions.html` explains the way to obtain the super regions suggested in `gauged_sites`. 

## Trend

The package `floodnetRfa` focuses on the analysis of station where flood risk is constant in time (stationary). 
This section performs standard trend tests.
The files `trend_script.R` evaluates the test statistics and the file `trend_readme.html` summaries them.



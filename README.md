# floodnetRfa extra

This repository includes supplemental material that comes along with the R-package `floodnetRfa`. It includes codes and summaries of the outcomes for the R-package when applied to the 1114 Canadian stations found in the dataset `gaugedSites` of `floodnetRfa`.
This material can be installed using the R terminal.

    library(devtools)
    install_github('floodnetProject16/floodnetRfa_extra')

The different folders are explained in the upcoming sections. 

## autopot

One method available to perform at-site flood frequency analysis is the method of Peaks Over Threshold (POT). 
The objective of this section is to use the automatic selection method provided in `floodnetRfa` to obtain candidates thresholds for all stations. 
The codes provided `autopot_script.R` serves to evaluate the thresholds provided in `gaugedSites`. A summary of the results is presented in the file `autopot_readme.html`.

## RFA

At-site and regional frequency analysis are performed on each of the 1114 stations. The files `atsite_script.R` and `RFA_script.R` fit the extreme distributions of the correspondings models.
Afterward, the file `comparison.html` compares the precision of the estimators based on the length of the confidence intervals (with coverage probability 95%). 

## Super regions

The recommended strategy for forming pooling groups first identifies a set of stations having similar properties called super regions. The file `super_regions.html` describes the way the super regions suggested in `gaugedSites` are obtained. 

## trend

The package `floodnetRfa` focuses on the analysis of station where constant flood risk is a reasonable hypothesis (stationary). This section performs common trend tests on the 1114 stations to verify the presence of trends.  The files `trend_script.R` evaluates the test statistics and the file `trend_readme.html reports a summary of the outcomes.

## Cache

Most of the previous scripts use parallel computing to speed up computation. To this end, a cache folder is created to save individual results for each site. For instance, `RFA/RFA_script.R` uses the `cache/RFA` folder. The detail of the files found in these folders is explained in the documentation found in the respective sections.
Comparison of the precison of flood frequency estimators
================
Martin Durocher
02/12/2019

# Introduction

The R-package `FloodnetRfa` provides functions to perform Regional
Frequency Analysis (RFA) following the guidelines of Floodnet. The
script `atsite_script.R` and `RFA_script.R` performed the analysis on
the 1114 stations using at-site and regional methods. This document
mainly compares the precision of the flood quantile estimates performed
by the different functions of `floodnetRfa`. The heterogeneity of the
pooling groups is also investigated.

# Evaluation criteria

The precision of two estimation methods is compared based on the length
of their confidence intervals (CIL), which is evaluated as the
difference between its upper and lower bounds on the interval with a
probability coverage of 95%. If \(r_i\) is the CIL of a reference method
and \(s_i\) is the CIL of an alternative method, we define the relative
difference of precision (RDP) as \((s_{i}-r_i)/r_i\). The average RDP
summarizes the overall difference in precision between two estimators
for a group of stations. Such differences in the RDP may be small and
simply due to the random nature of the experiment. To assess the
significance of the difference in precision, we evaluate the p-value of
the paired Wilcoxon signed-rank test (WILCOX). This test assumes the
comparison of symmetric distributions around the median, which was found
to be a reasonable hypothesis after a logarithm transformation. In
particular, a p-value below the significance level of `0.05` will be
considered as an indication of a significant different of precision.

Another way of comparing two estimation methods consists to examine the
percentage of stations (PCT) where the CIL of the alternative estimator
is lower than the one of the reference estimator (\(s_i < r_i\)). In
this case, the p-value of a binomial test (PVAL) allows verifying the
hypothesis that the PCT is significantly different from 50%. One
drawback of this approach is that the magnitude of the difference is not
considered.

# At-site Flood Frequency Analysis

## FloodnetAmax vs FloodnetPot

In a first step, we investigated at-site frequency analysis for both
Annual Maximum discharge (AMAX) and Peaks Over Threshold (POT). In other
words, we compares the precision of the function `FloodnetAmax` versus
`FloodnetPot`, where the treshold of the POT model is automatically
selected. The results are produced by the script `atsite_script.R`,
which produces summary graphics that are saved in the CACHE folder. An
example is presented at the end of the document.

First, we noticed that in some situations, POT performed poorly.
Therefore, we examine the stations where the CIL of the POT estimates is
more than 10 times the one of AMAX. These situations are consider as
clear cases where AMAX outperforms POT. Such situations occur for 16
stations, of which 13 are found in region 3.

Table 1 reports the evaluation criteria described earlier. Both the
WILCOX and Binomial tests suggest a significant difference in favour of
POT for Q100. One exception is regions 3 where WILCOX indicates that the
difference is not significant. Regions 3 have a positive RDP, but a PCT
higher than 50% for POT, which indicates that POT is more precise in
most situations, but tends to be substantially less precise in the
others. According to the binomial test, the difference in precision for
Q10 is not significant for regions 1, 4, 5 and 6. On the opposite, POT
is more often more precise in region 3.

It should be considered that for region 3, several stations are seasonal
stations that don’t record streamflow all year long. To preform the
analysis, the missing days were assumed as values below the threshold.
For this practical reason, AMAX should be recommended as the estimation
is based on a strong assumption. Overall, the results suggest no clear
evidence for Q10 that one method is generally more precise than the
other. On the other hand, POT appears more often more precise when
considering Q100.

**Table 1: Comparison of precision between FloodnetAmax and
FloodnetPot.**

***Q10***

    ##           1    2    3    4    5    6  tot
    ## RDP    -0.6 -6.7 -4.6  0.9  5.0 -2.1 -1.2
    ## WILCOX  6.4  0.0  0.0  2.8 80.9  1.2  0.0
    ## PCT    58.6 68.0 63.9 56.0 52.0 55.0 58.6
    ## PVAL    7.7  0.0  0.0  6.6 65.0 20.5  0.0

***Q100***

    ##            1    2    3    4    5     6  tot
    ## RDP    -22.0 -7.6 34.3 -4.8  2.3 -19.0  1.9
    ## WILCOX   0.0  0.0 45.1  0.0  0.0   0.0  0.0
    ## PCT     81.9 63.9 57.1 70.4 62.9  79.4 67.9
    ## PVAL     0.0  0.8  2.0  0.0  0.1   0.0  0.0

# Regional frequency analysis

## FloodnetAmax vs FloodnetPool

In this section, we compare the estimation of the at-site
(`FloodnetAmax`) and regional method (`FloodnetPool`) using the
L-moments algorithm and annual maximum streamflow. The Regional
Frequency Analysis (RFA) is using a pooling group approach starting with
a region of influence including 25 stations that is later updated by
removing the most heterogenous stations.

Using the 1114 stations, Table 2 compares the CIL of the at-site and
regional methods. The criteria show that RFA improves the precision of
the estimates for a large proportion of the stations. For Q100, RFA
leads to more precise estimates in 87.5% of the stations. The region
with the lowest improvement is region 4, with an RDP of -6.9%, while
across Canada an RDP of -30.9% is observed. Note that all comparisons
are significants.

**Table 2: Comparison of precision between FloodnetAmax and
FloodnetPool**

***Q10***

    ##            1     2    3    4     5     6  tot
    ## RDP    -12.8 -14.3 -1.6 -2.7 -11.9 -15.5 -8.0
    ## WILCOX   0.0   0.0  0.0  0.0   0.0   0.0  0.0
    ## PCT     87.9  83.5 62.8 70.0  86.9  85.6 76.3
    ## PVAL     0.0   0.0  0.0  0.0   0.0   0.0  0.0

***Q100***

    ##            1     2     3    4     5     6   tot
    ## RDP    -38.9 -30.8 -30.6 -6.9 -45.1 -46.0 -30.9
    ## WILCOX   0.0   0.0   0.0  0.0   0.0   0.0   0.0
    ## PCT     92.2  84.5  89.8 73.9  94.9  94.4  87.5
    ## PVAL     0.0   0.0   0.0  0.0   0.0   0.0   0.0

## FloodnetPool: AMAX vs POT

Next, Table 3 compares the precision RFA-AMAX and RFA-POT methods that
are both performed by `FloodnetPool`. For regions 1 to 3, the binomial
and WILCOX test suggest that the precision of Q100 is not significantly
different between the two methods. For regions 5 and 6, AMAX is
significantly more precise with RDP of 21.9% and 17.9%, while POT is
more precise in 67.6% of the stations in region 4.  
Beside region 3, RFA-AMAX has better RDP and PCT for all regions.

**Table 3. Comparison of precision FloodnetPool using AMAX and POT**

***Q10***

    ##           1    2    3    4    5    6  tot
    ## RDP    20.4  5.6 -2.6  5.8 14.6 19.9  8.8
    ## WILCOX  0.0  8.7  0.0  7.8  0.0  0.0  0.0
    ## PCT    25.9 45.4 62.8 47.0 22.9 17.8 40.3
    ## PVAL    0.0 41.7  0.0 37.9  0.0  0.0  0.0

***Q100***

    ##           1    2    3     4    5    6  tot
    ## RDP     1.8 -3.4  4.1 -10.4 21.9 17.9  4.9
    ## WILCOX 12.4  6.9  5.2   0.0  0.0  0.0 19.8
    ## PCT    57.8 53.6 53.2  67.6 20.0 30.0 48.0
    ## PVAL   11.4 54.3 29.3   0.0  0.0  0.0 19.8

# Overall scores

Finally, we summarize the relative precision of all methods in terms of
a single score. To do so, we rank the CIL from the least to the most
precise by estimation method. Next, we assign a score from 0 to 3 and
evaluate the average score by regions.

Table 5 indicates the sum of the scores for Q10 and Q100 as an overall
assessment of the modeling precision. The resulting statistics can take
values between 0 and 6. As expected from the pairwise comparison, the
at-site methods have lower scores than the regional methods and the
at-site POT method does better than the at-site AMAX method. Overall,
the regional method using AMAX has the best scores for regions 1, 2, 5
and 6, while the regional POT methods has the best scores for regions 3
and 4.

**Table 5: Comparison of all functions based on precision scores.**

    ##            1   2   3   4   5   6 all
    ## at-amax  1.3 1.0 1.5 1.4 1.4 1.1 1.3
    ## at-pot   2.8 2.3 2.3 2.7 1.8 2.2 2.4
    ## rfa-amax 4.4 4.5 3.6 3.8 5.1 5.1 4.4
    ## rfa-pot  3.5 4.2 4.5 4.1 3.8 3.5 3.9

# Heterogeneity

``` r
het <- read.csv(gzfile('heteo.csv.gz'))

CutH <- function(z){
    cut(z, c(-Inf,1.5,2, Inf), 
            labels = c('Homo.','Prb.','Hetero.'))
}

CutN <- function(z){
    cut(z, c(-Inf,30,60, Inf), 
            labels = c('Short','Medium','Long'))
}

hamax <- CutH(het$hamax)
namax <- CutN(het$namax)
hpot <- CutH(het$hamax)
npot <- CutN(het$namax)
```

**Table {\#fig:hetero} : Heterogeneity of the pooling groups by super
regions.**

test @ref(fig:hetero)

``` r
round(100*addmargins(table(H_AMAX = hamax, Region = supreg))/1114,1)
```

    ##          Region
    ## H_AMAX        1     2     3     4     5     6   Sum
    ##   Homo.     2.4   1.0   6.6   2.8   6.2   8.0  26.9
    ##   Prb.      5.7   6.5  13.4  10.4   9.5   8.0  53.5
    ##   Hetero.   2.2   1.3   6.4   9.5   0.0   0.2  19.6
    ##   Sum      10.4   8.7  26.3  22.7  15.7  16.2 100.0

``` r
round(100*addmargins(table(H_POT = hpot, Region = supreg))/1114,1)
```

    ##          Region
    ## H_POT         1     2     3     4     5     6   Sum
    ##   Homo.     2.4   1.0   6.6   2.8   6.2   8.0  26.9
    ##   Prb.      5.7   6.5  13.4  10.4   9.5   8.0  53.5
    ##   Hetero.   2.2   1.3   6.4   9.5   0.0   0.2  19.6
    ##   Sum      10.4   8.7  26.3  22.7  15.7  16.2 100.0

# Appendix

## Summary graph of the at-site analysis.

The graphic below is a summary produced by `at-site.R`. The left panels
show the return level plots of the at-site AMAX and POT method for
stations `01AD002`. The number of observations (N) is indicated below
the x-axis. The top right panel presents the estimated flood quantiles
with the confidence interval (95%) and the bottom right panel shows the
respective coefficient of variations, *i.e.* the standard deviation
standardized by the flood quantile.

<center>

<img src="atsite_01AD002.png" height="960" width="900">

</center>

## Summary Graphic of the regional analysis.

The graphic below is a summary produced by `RFA.R`. The top panels
illustrate the pooling group of a target station (in red) for the
seasonal and descriptor space. The members of the pooling group are in
blue and the black dots are the other stations. The bottom panels
present the estimated flood quantiles with the confidence interval (95%)
by return periods as well as the coefficient of variations.

<center>

<img src="rfa_01AD002.png" height="960" width="900">

</center>

## Super regions

Results are reported by regions having similar geographical locations
and site characteristics to help to generalize conclusions. Figure 1
illustrates the descriptor space defined by the drainage area (AREA) and
Mean Annual Precipitation (MAP). Here is a brief description of the
considered regions.

  - Super region 1 are stations found mainly in the northern part of
    Canada. The basins are generally large with few annual
    precipitations.

  - Super region 2 are stations found mainly on the Pacific coast and
    the Fraser Valley.  
    The basins have important annual precipitation.

  - Super region 3 are stations found mainly in the Prairies and the
    mountains of British Columbia. The basins are generally small with
    few annual precipitations.

  - Super region 4 are stations found mainly in the north of the
    Canadian western provinces. The basins are generally large with
    important annual precipitations.

  - Super region 5 are stations found mainly in Southern Ontario. The
    basins are generally small with important annual precipitations.

  - Super region 6 are stations found mainly in the eastern part of
    Canada. The basins generally have important annual precipitations.

![](RFA_readme_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

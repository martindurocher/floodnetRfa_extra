# Meta-information about HYDAT gauged stations.

**_File: gauged_sites.csv_**

The dataset contains meta-information about 1114 stations in HYDAT that have at least 20 years of observations and a natural flood regime.
This information can be used to perform flood frequency analysis using either the annual maximum (AMAX) or the peaks over threshold (POT) methods. 
The dataset contains basic information about the stations and widely available catchment descriptors.

* Identification number: `station`. 
* Province: `prov`.
* Coordinates:  `lon`, `lat`.
* Drainage area: `area`.
* Mean Annual Precipitation `map`.

A super region is a group of stations that have similar properties and serves as an initial population to form a pooling group based on seasonal statistics.
The dataset suggests four groups of super regions that come from the application of the Ward's (`hc`) and k-means (`km`) clustering techniques.
For example, the column `supreg_km12` is the result of the classification of the stations into 12 super regions by the k-means algorithm.
The seasonal statistics characterize the regularity (`season_radius`) and timing  (`season_angle`) of the annual floods. 
These values represent polar coordinates inside the unitary circle and are also available as cartesian coordinates (`season_x`, `season_y`).

The dataset includes candidates thresholds for POT analysis for each station.
The candidates in the column `auto` correspond to the first thresholds that have a p-value greater or equal to 0.25 and an exceedance rate lower than 2.5 peaks per year (PPY).
The other suggested thresholds are associated with given PPY.
For example, `pp175` corresponds to approximately 1.75 PPY.

The dataset contains p-values of standard trend tests to identify stations with time-dependent flood risk. 
For AMAX, it includes the outcomes of the Mann-Kendall's (`trend_mk`) and the Pettitt's test (`trend_pt`).
For POT, it includes the outcomes of the Mann-Kendall's test applied to the exceedances (`trend_mx`) and of logistic regression models (`trend_lg`). 
More precisely, F-tests compare polynomial trends of order 1 to 3 to a constant logistic regression model of the exceedance rate. 
The lowest p-value serves as evidence that at least one of the considered trends was significant.

# Catchment descriptors

**_File: descriptors.csv_**

The dataset contains a list of meteorological and physical characteristics of 770  gauged stations. This dataset is useful for estimating flood quantiles at ungauged sites. It includes the following variables.

* Station identification number: `station`.
* Station geographical coordinates: `lon`, `lat`.
* Watershed centers: `lon_ws, lat_ws`.
* Drainage area (m3/s): `area`.
* Perimeter (km): `peri`.
* Station elevation (m): `elev`.
* Watershed average elevation (m): `elev_ws`.
* Mean aspect (Deg): `aspect`.
* Mean slope (%): `slope`.
* Stream density ($km^{-1}$): `stream`.
* Percentage of water bodies (%): `wb`.
* Station and watershed mean annual precipitation (mm): `map`, `map_ws`.
* Station and watershed annual average air temperature (C): `temp`, `temp_ws`.

# References

1. Durocher, M., Burn, D. H., & Mostofi Zadeh, S. (2018). A nationwide regional flood frequency analysis at ungauged sites using ROI/GLS with copulas and super regions. Journal of Hydrology, 567, 191-202. https://doi.org/10.1016/j.jhydrol.2018.10.011

2. Durocher, M., Zadeh, S. M., Burn, D. H., & Ashkar, F. (2018). Comparison of automatic procedures for selecting flood peaks over threshold based on  goodness-of-fit tests. Hydrological Processes, 0(0). https://doi.org/10.1002/hyp.13223

3. Mostofi Zadeh, S., & Burn, D. H. (2019). A Super Region Approach to Improve Pooled Flood Frequency Analysis. Canadian Water Resources Journal / Revue Canadienne Des Ressources Hydriques, 0(0), 1–14. https://doi.org/10.1080/07011784.2018.1548946
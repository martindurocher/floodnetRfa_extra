##################################################################
# I used classical multidimensional scaling (CMD) on some occasions to transform 
# coordinates into a Cartesian system that preserves the great-circle distance. 
# This code verifies that this approach compares well with  
# other alternatives, such as Lambert's conformal conic projection (LCC) and 
# non-metric multidimensional scaling (NMD). 
##################################################################

library(floodnetRfa)
library(CSHShydRology)
library(sp)
library(MASS)

source('config')

## Compute the true Great-Circle distances
coord <- GAUGEDSITES[, c('lon','lat')] 
h <- GeoDist(coord)
v <- h[lower.tri(h)] 

## Compute the LCC
xdesc <- GAUGEDSITES[, c('lon','lat','station')]
coordinates(xdesc) <- ~lon+lat
proj4string(xdesc) <- CRS("+init=epsg:4326")
xdesc <- spTransform(xdesc, CRS("+init=epsg:3347"))

## Evaluate distances
hlcc <- as.matrix(dist(coordinates(xdesc)))/1000
vlcc <- hlcc[lower.tri(hlcc)]

## Compute the CMD and respective distances
coord2 <- -cmdscale(h)
hcmd <- as.matrix(dist(coord2))
vcmd <- hcmd[lower.tri(hcmd)]

## Compute the
## add a small pertubation to very close sites. Avoid distance = 0. 
h[426,425] <- h[425,426] <- .001
coord3 <- isoMDS(h)
hnmd <- as.matrix(dist(coord3$points))
vnmd <- hnmd[lower.tri(hnmd)]

## Compate the results using absolute differences
summary(abs(vlcc - v))
summary(abs(vcmd - v))
summary(abs(vnmd - v))

dif.lcc <- abs(vlcc - v) - abs(vcmd - v)
wilcox.test(dif.lcc)
dif.nmd <- abs(vnmd - v) - abs(vcmd - v)
wilcox.test(dif.nmd)

#################################################################
## Conclusion
##
## On one side, the difference between CMD and NMD is marginal
## and CMD is faster.
## On the other hand CMD is shown to better preserve the distances
## than using the LCC. 
#################################################################
 
############################################################
## Objective is to save a matrix of the intersite correlation
## coefficient for all sites in gauged_sites.csv that could be used later on.
## The resulting matrix must be postive definite (PD) if we want to use it for 
## simulation. In particular, we must be able to perform the Cholesky 
## decomposition that is often used for simulating multivariate normal samples.
############################################################

library(floodnetRfa)
library(CSHShydRology)
library(sp)

rm(list = ls())

## Path of the HYDAT database
## Create a variable DB_HYDAT
source('config')

## utility function
ToVec <- function(z) z[lower.tri(z)]

## Compute the true Great-Circle distances
coord <- GAUGEDSITES[, c('lon','lat')] 
h <- GeoDist(coord)
vec.h <- ToVec(h) 
vec.near <- (vec.h < 500)

## Read data and format them in a matrix with sites in columns
xall <- AmaxData(DB_HYDAT, GAUGEDSITES[,'station'])
xall$year <- as.integer(format(xall$date,"%Y"))
xw <- DataWide(value ~ site + year, xall)

## Compute the number of paired observations
nmat <- crossprod(!is.na(xw))

##################################################################
## Compute the the Spearman rhos
cc.raw <- cor(xw, method = 'spearman', use = "pairwise.complete.obs") 

## Remove paired observation with too few pairs
cc.raw[nmat < 10] <- NA

## Convert from Spearman to Pearson coefficient
cc.raw <- 2*sin(pi/6*cc.raw)

###################################################################
# Correct correlation matrix to be PD.
cc.emp <- Intersite(xw)$model

## verify that it can be decomposed
o <- chol(corr_emp)

##############################################################
## Using exponential model
cc.exp <- Intersite(xw, type = 'exp', distance = h, smooth = .6)$model

## verify that it can be decomposed
o <- chol(corr_emp)

###############################################################
## Visualized the results
###############################################################

# Limit data to a smaller sample for visualization
# and extract paired observations
sid <- sample(which(vec.near), 2000)
h0 <- vec.h[sid]

ToVec0 <- function(z) ToVec(z)[sid]
pair.raw  <- ToVec0(cc.raw)
pair.emp  <- ToVec0(cc.emp)
pair.exp  <- ToVec0(cc.exp)

## Fitting of the exponential model on raw coefficient
plot(h0, pair.raw)
points(h0, pair.exp, col = 'orange', pch = 16, cex = .5)
abline(h=0, col = 'cyan', lwd = 3)

## Same with PD corrected coefficient
plot(h0,pair.emp)
points(h0, pair.exp, col = 'orange', pch = 16, cex = .5)
abline(h=0, col = 'cyan', lwd = 3)

## Difference between raw and corrected coefficient
dif <- pair.emp - pair.raw
dd <- !is.na(dif)

plot(h0,dif); abline( h =0, col = 'cyan', lwd = 3)
lines(smooth.spline(h0[dd],dif[dd]), col = 'orange', lwd = 3)

## Descriptive statistics for discrepencies
summary(dif[dd]) ## avg = -0.03, IQR = 0.086
sd(dif[dd])  ## 0.083

dif <- pair.exp - pair.raw
dd <- !is.na(dif)
summary(dif[dd]) ## avg = -0.03, IQR = 0.086
sd(dif[dd])  ## 0.083

###############################################################################
## NOTE:
## Overall the corrected coefficients preserve the correlation decay as a function 
## of the distance. We see that the PD correction has a small tendency to 
## systematically lower the correlation coefficients in comparison with the raw 
## estimates. The effect is more pronounced for closer sites. The 
## differences are reasonably low in most cases, with an IQR of 0.08. In 
## comparison, the exponential model has an IQR of 0.30.
###############################################################################

## Uncomment to save new results
outfile <- gzfile('./other/intersite_emp.csv.gz')
#write.csv(corr_emp, file = outfile, row.names = FALSE)

outfile <- gzfile('./other/intersite_exp.csv.gz')
#write.csv(corr_exp, file = outfile, row.names = FALSE)

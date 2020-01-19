library(floodnetRfa)
library(foreach)
library(doParallel)

rm(list = ls())

## Load function for the analysis
source('trend/trend_lib.R')

## Path of the HYDAT database
## Create a variable DB_HYDAT
source('config')

## Size of the bootstrap sample
NSIM <- 2000

## Resolution of the output figure
RES <- 560

## Polynomial degree consider for smoothing graph should be 1:3
DEG <- 1:3

## Number of parallel processus
NCPU <- 4

## Temporary files used during parallel computing
DIR_CACHE <- 'cache/trend'

################################################
## BODY
################################################

if(!file.exists('cache'))
	dir.create('cache')

if(!file.exists(DIR_CACHE))
	dir.create(DIR_CACHE)

progress.file <- file.path(CACHE, 'progress.log')

## Load the list of stations to analyse
stations <- gaugedSites$station

## Monitoring file for parallel computing
t0 <- Sys.time()
write(format(t0), file = progress.file)

## Set-up parallel computing
cl <- makeCluster(NCPU)
registerDoParallel(cl)

statu <- foreach(ii = seq_along(stations),
                 .errorhandling = 'pass',
                 .packages = c('CSHShydRology', 'floodnetProject16',
                 							'Kendall','boot','trend')) %dopar%{

  site <- stations[ii]

  write(paste0('(', ii ,') Analyzing Trend ', site),
	      file = progress.file, append = TRUE)

  ##########################
	## Trend in the AMAX data
	##########################

  ## read annual maximum from HYDAT
  an <- AmaxData(site,DB_HYDAT)[,-1]
  nyear <- nrow(an)

  ## Compute the p-value of the MannKendall test
  mk <- AutoMK(an$value)
  mk.pvalue <- mk[2]

  ## Compute the p-value of the Pettitt test
  pt <- pettitt.test(an$value)
  pt.pvalue <- pt$p.value

  ##########################
	## Trend in the POT
	##########################

  ## Read Daily data
  xd <- DailyData(site, DB_HYDAT)[,-1]

  ## Identify peaks using previously selected thresholds
  u <- gaugedSites[gaugedSites$station == site, 'auto']
  area <- gaugedSites[gaugedSites$station == site, 'area']
  peaks <- which.floodPeaks(value~date, xd , u = u, r = 4 + log(area))

  ## Compute p-value for the logistic model (Trend in number of peaks)
  logis.pvalue <- TrendLogis(xd$date, peaks)

  ## Compute p-value for the Mann-Kendall test (Trend in mean excess)
  mex.pvalue <- AutoMK(xd[peaks,'value'])[2]

  ####################################
  ## save the p-values
  ####################################
  all.pvalue <-
  	data.frame(trend_mk = mk.pvalue,
  						 trend_pt = pt.pvalue,
  						 trend_lg = logis.pvalue[1],
  						 trend_mx = mex.pvalue)

  pval.file <- paste0(paste0(site, '.csv'))
  write.csv(round(all.pvalue,4),
  					file = file.path(CACHE, pval.file),
  					row.names = FALSE)

  ####################################
  ## Summary graphics
  #####################################

  ## Compute the fitted value of the change point identify in the Pettitt test
  ptk <- pt$estimate
  pt1 <- mean(an$value[1:ptk])
  pt2 <- mean(an$value[(ptk +1):nyear], na.rm = TRUE)
  pt.yhat <- rep(pt2,nyear)
  pt.yhat[1:ptk] <- pt1


  fig.file <- paste0(paste0(site, '.png'))
	png(file = file.path(CACHE, fig.file), height = 1.33*RES, width = RES)

	par(mfrow = c(2,1), mar = c(5,5,5,5))

	##---- Plot the AMAX ----##
  plot(an, type = 'b',
  		 xlab = 'Year', ylab = 'Annual maxima',
  		 main = paste0(site,': MK=', round(mk.pvalue,3),
  		 							 '; PT=', round(pt.pvalue,3) ))

  ## estimate a smooth trend
  FunFit <- function(z) lm(log(value) ~ poly(year,z), an)
  an.sm <- Map(FunFit, DEG)
  an.sm <- an.sm[[which.min(sapply(an.sm,AIC))]]

  lines(an$year, exp(predict(an.sm)), col = 'red', lwd = 3)

  ##----- Plot of the POT -----##
  xd.peaks <- xd[peaks,]
  plot(xd.peaks, type = 'h', ylab = 'Flood peaks', xlab = 'Date',
  		 main = paste0('Trend : MX=', round(mex.pvalue,3),
  		 							 '; LG=', round(logis.pvalue,3) ))

  ## Add a smooth lines for the mean excess
  FunFit <- function(z) lm(log(value)~poly(date,z), xd.peaks)
  xd.sm <- Map(FunFit, DEG)
  xd.sm <- xd.sm[[which.min(sapply(xd.sm,AIC))]]
  lines(xd.peaks$date, exp(predict(xd.sm)), col = 'red', lwd = 3)

  ## Add a smooth line for the number of PPY
  npy <- tapply(seq_along(xd$date) %in% peaks, format(xd$date,'%Y'), sum)
  npy <- data.frame(value = npy,
  									date = as.Date(paste0(names(npy),'/7/15')))

  FunFit <- function(z) glm(value~poly(date,z), npy, family = poisson())
  np.sm <- Map(FunFit, DEG)
  np.sm <- np.sm[[which.min(sapply(np.sm,AIC))]]

  par(new = TRUE)
  plot(npy$date, predict(np.sm, type = 'response'),
  		 type = 'l', axes = FALSE, ylab = NA, xlab = NA,
  		 pch = 16, col = 'blue', ylim = c(0,5), lty = 2, lwd = 3)

  axis(side = 4)
  mtext(side = 4, line = 3, 'Number of PPY')

  legend('top', horiz = TRUE, legend = c('Mean excess','trend in PPY'),
  			 lty = 1:2, col = c('red','blue'))

  dev.off()

  return(0)
}

stopCluster(cl)

###########################################
## Read and merging the output
###########################################

sfiles <- list.files(CACHE, full.names = FALSE, pattern = '*.csv')
snames <- substr(sfiles,1,7)
trend <- lapply(file.path(CACHE, sfiles), read.csv)
trend <- cbind(station = snames, do.call(rbind, trend))

## verify which stations did not succeeded
fail.id <- which(!gaugedSites$station %in% trend$station)
print(fail.id)

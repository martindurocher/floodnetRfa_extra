library(floodnetRfa)
library(foreach)
library(doParallel)

rm(list = ls())

## Path of the HYDAT database
## Create a variable DB_HYDAT
source('config')

## Resolution of the graphics in pixels (RES x RES)
RES_FIG <- 640

## Number of parallel processes
NCPU <- 4

## Cache folder
DIR_CACHE <- 'cache/autopot'

##########################################
## BODY
##########################################

if(!file.exists('cache'))
	dir.create('cache')

if(!file.exists(DIR_CACHE))
	dir.create(DIR_CACHE)

progress.file <- file.path(DIR_CACHE,'progress.log')

## Load the list of stations to analyse
stations <- gaugedSites[,1]

## Determine the minimal time between extracted peaks
rarea <- ceiling(4 + log(gaugedSites[,'area']))

## Set-up parallel computing
cl <- makeCluster(NCPU)
registerDoParallel(cl)

## Monitoring progress for parallel computing
t0 <- Sys.time()
write(format(t0), file = progress.file)

statu <- foreach(ii = 1:8, #seq_along(stations),
								 .packages = c('CSHShydRology','floodnetRfa'),
								 .errorhandling = 'pass') %dopar%{

	site <- stations[ii]

	write(paste0(ii,' - Analyzing ', site), file = progress.file, append = TRUE)

  ###########################################
	## Search for thresholds
	###########################################

	xd <- try(DailyData(site, DB_HYDAT, pad = TRUE, tol = 346)[,-1])

	nyear <- ifelse(class(xd) == 'try-error', 0,
				          length(unique(format(xd$date, '%Y'))))

	## In case of seasonal station
	if((nyear < 15)){
		xd <- DailyData(site, DB_HYDAT, pad = TRUE, tol = 60)[,-1]
		nyear <- length(unique(format(xd$date, '%Y')))
	}

	## pre-select a set of candidate thresholds based on observations
	candidate <- sort(unique(xd$value), decreasing = TRUE)

	## allocate memory
	umat <- matrix(0, length(candidate), 6)
	colnames(umat) <- c('u','PPY','AD','Shape', 'MRL', 'Q10')

	## Compute the p-value of the Anderson Darling test for all candidates
	## plus other statistics
	for(jj in seq_along(candidate)){

		suppressWarnings(fit0 <- try(
			FitPot(value~date, xd, u = candidate[jj], declust = 'wrc', r = rarea[ii]),
		  silent = TRUE))

		## If POT fails
    if(class(fit0) != 'fpot'){
      umat[jj,] <- c(candidate[jj], rep(0,ncol(umat)-1))

    ## If less than 20 peaks
    } else if(fit0$nexcess < 20){
      umat[jj,] <- c(candidate[jj], rep(0,ncol(umat)-1))

    ## Normal condition perform the AD test and compute the flood quantile Q10
    } else{

    	hat0 <- predict(fit0, rt = 10)
      gof0 <- GofTest(fit0)$pvalue

      umat[jj,] <- c(candidate[jj],
      							 fit0$nexcess/fit0$nyear,
      							 gof0,
      							 fit0$estimate[2],
      							 fit0$mrl,
      							 hat0)
    }

		## Terminate loops when PPY is greater than 4
    if(umat[jj, 'PPY'] > 4)
    	break

	}#endfor

	## remove fails
	umat <- as.data.frame(umat)
	umat <- umat[umat$PPY > 0 & umat$u > 0, ]

	##########################################
	## Find thresholds of given PPY
	##########################################

	## Because of the declustering techniques the association between
	## PPY and candidate thresholds (u) is not a strictly monotone
	## We look for a smooth approximation to obtain the threshold associated
	## to a given PPY.

	## PPY of interest: 1, 1.25, ..., 2.5
	ppy.step <- seq(1, 2.5, 0.25)

  if(nrow(umat) < 1){
  	## Case there is no value to evaluate
  	## use the minimum value
	  uref <- rep(min(xd$value), length(ppy.step))

	} else {

		## isotonic regression is used to smooth the PPY
		xppy <- try(data.frame(u = umat$u, PPY = isoreg(umat$PPY)$yf))

		ppy.adiff <- abs(sapply(ppy.step,'-', xppy$PPY))
		ppy.id <- apply(ppy.adiff, 2, which.min)
		uref <- rev(xppy$u[ppy.id])

	}

	names(uref) <- paste0('ppy',rev(ppy.step)*100)

	#############################################
	## Identify best threshold
	#############################################

	## Identify an automatic threshold.
	if(nrow(umat) < 1){
	  u0 <- min(xd$value)
	} else {

	  u25 <- umat[umat$u >= uref[1], ]
	  usgn <- u25[u25$AD >=.25, ]$u
	  umax <- u25[u25$AD >= max(u25$AD), ]$u

	  u0 <- ifelse(length(usgn) > 0, min(usgn), min(umax))

	}

	#############################################
  ## Seasonal analysis
  #############################################

	an <- try(AmaxData(site, DB_HYDAT, year = FALSE), silent = TRUE)

	if(class(an) == 'try-error'){
	  an <- ExtractAmax(value~date, xd, tol = 60)
	}

	## For annual maxima
	if(nrow(an) > 5){
  	ss <- SeasonStat(an$date)
	} else {
		ss <- c(x = 0, y = 0, angle = 0, radius = 0)
	}

	names(ss) <- paste0('season_', names(ss))

  ######################################
	## graphics for validating the threshold
	######################################

	if(nrow(umat) > 1){

		## Get the range of value for plot limits
  	rg <- apply(umat, 2, range)
	  rg.diff <- apply(rg,2,diff) *.1

    fname <- file.path(DIR_CACHE, paste0(site,'.png'))
	  png(fname,  width = 3*RES_FIG, height = 2*RES_FIG, pointsize = 24)

	  palette(c('#000000','#e41a1c','#4daf4a','#377eb8','#984ea3',
		  			'#ff7f00','#a65628','#f781bf'))

	  even.id <- seq(1,length(uref),2)

	  par(mfrow = c(2,2))

	  ## Plot main label
    mlabs <- c('','',
  					   paste0('Station ID = ',site),
  					   paste0('Auto threshold, u = ', round(u0, 1)),
  				  	 paste0('Nb. year = ', round(nyear,1)),
  			  		 paste0('Season(month,radius) = (', round(ss[3]/pi*6) ,', ',
  		  			 			 round(ss[4],2), ')')
  	  				 )

	  for(kk in 3:6){
	    plot(umat[,c(1,kk)], pch = 16, type = 'l',
	  		   ylim = rg[,kk] + rg.diff[kk] * c(-1,1),
	  		   main = mlabs[kk],
	  	  	 xlab = 'Threshold (u)',
	    		 ylab = colnames(umat)[kk], log = 'x')
	    abline(v = u0, col = 2, lwd = 3)
	    abline(v = uref, lty = 3, col = 5)

	    text(x = rev(uref)[even.id], y = rg[1,kk] - rg.diff[kk]/2,
	  	  	 labels = ppy.step[even.id],
	    		 pos = 1, cex = 1.2, col = 5)
	  }

	  dev.off()
	}
	####################################
	## Output the threshold values
	####################################

	fname <- file.path(DIR_CACHE, paste0(site,'.csv'))
	out <- list(paste0('"',site,'"'),
								nrow(an),
								round(nyear,1),
								round(ss,3),
								round(u0,1),
								round(uref,1))

	write(paste(unlist(out), collapse = ','), file = fname)

	return("exit normally")

}#endforeach

stopCluster(cl)

write(format(t0-Sys.time()),
			file = progress.file, append = TRUE)

statu <- lapply(statu,'[',1)
print(table(unlist(statu)))

###########################################
## Read and merging the output
###########################################

lst <- list.files(DIR_CACHE, '*.csv')

Fz <- function(z) read.csv(file = file.path(DIR_CACHE, z), header = FALSE)
xthresh <- lapply(lst,Fz)
xthresh <- do.call(rbind, xthresh)
colnames(xthresh) <- c('station', 'year.amax', 'year.pot',
											 'season_x', 'season_y', 'season_angle', 'season_radius',
											 'auto', 'ppy250', 'ppy225', 'ppy200',
											 'ppy175', 'ppy150', 'ppy125', 'ppy100')

## verify which stations did not succeeded
fail.id <- which(!gaugedSites$station %in% xthresh$station)
print(fail.id)

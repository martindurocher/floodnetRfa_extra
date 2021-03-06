require(floodnetRfa)
require(foreach)
require(doParallel)

rm(list = ls())

## Path of the HYDAT database
## Create a variable DB_HYDAT
source('config')

## Size of the bootstrap sample
NSIM <- 2000

## Number of parallel processes
NCPU <- 3

## Resolution figures
RES <- 480

## script and cache folder
DIR_OUT <- 'RFA'
DIR_CACHE <- 'cache/atsite'

################################################################################
## Monitoring file for parallel computing

if(!file.exists('cache'))
	dir.create('cache')

if(!file.exists(DIR_CACHE))
	dir.create(DIR_CACHE)

progress.file <- file.path(DIR_CACHE,"progress.log")
error.file <- file.path(DIR_CACHE,"error.log")

t0 <- Sys.time()
write(format(t0), file = progress.file)
write(format(t0), file = error.file)

## Set-up parallel computing
cl <- makeCluster(NCPU)
registerDoParallel(cl)

statu <- foreach(ii = 1:4, #nrow(GAUGEDSITES),
								 .packages = 'floodnetRfa',
								 .errorhandling = 'pass') %dopar%{

	site <- as.character(GAUGEDSITES$station[ii])
	thresh <- GAUGEDSITES$auto[ii]
	area <- GAUGEDSITES$area[ii]

	####################
	## Carry out AMAX ##
	####################

  write(paste0(ii, ' - Analyzing ', site),
				file = progress.file, append = TRUE)

	an <- AmaxData(DB_HYDAT, site)
	amax <- try(FloodnetAmax(an, out.model = TRUE,
													verbose = FALSE, nsim = NSIM), silent = TRUE)

	if(class(amax) != 'try-error'){
		amax.file <- paste0('amax',site,'.csv')
		write.csv(amax$qua, file = file.path(DIR_CACHE, amax.file), row.names = FALSE)
		out.amax <- TRUE

	} else {
		write(paste0('Error analyzing AMAX for station: ', site,']'), 
					file = error.file, append = TRUE)
		out.amax <- FALSE
	}

	###################
	## Carry out POT ##
	###################

	xd <- try(DailyData(DB_HYDAT, site, pad = TRUE, tol = 346)[,-1])

	nyear <- ifelse(class(xd) == 'try-error', 0,
				          length(unique(format(xd$date, '%Y'))))

	## In case of seasonal stations
	if((nyear < 15)){
		xd <- DailyData(site, DB_HYDAT, pad = TRUE, tol = 60)[,-1]
		nyear <- length(unique(format(xd$date, '%Y')))
	}


	pot <- try(FloodnetPot(xd, u = thresh, area = area,
												 verbose = FALSE, nsim = NSIM,
												 out.model = TRUE), silent = TRUE)

	if(class(pot) != 'try-error'){
		pot.file <- paste0('pot',site,'.csv')
		write.csv(pot$qua, file = file.path(DIR_CACHE, pot.file), row.names = FALSE)
		out.pot <-  TRUE

	} else {
		write(paste0('Error analyzing POT for station: ', site,']'), 
					file = error.file, append = TRUE)
		out.pot <- FALSE
	}

	###################
	## Summary GRAPH ##
	###################


	fig.file <- paste0(paste0(site, '.png'))
	png(file = file.path(DIR_CACHE, fig.file), height = 1.5*RES, width = 2*RES)

	layout(matrix(1:4,2,2))

	## allocate memory
	rg.amax <- rg.pot <- NULL

	if(out.amax)
		plot(amax) + 
			labs(title = paste0('AMAX: ', site), caption = paste0('N: ',nb.amax)) 
		

	if(out.pot)
		plot(pot) + 
			labs(title = paste0('POT'), caption = paste0('N: ',nb.pot)) 

	cp <- CompareModels(amax, pot)
	
	plot(cp) + theme(legend.position = 'top') 
	plot(cp, 'cv') + theme(legend.position = 'top') 
	
	## Return level plot ##
	rg.qua <- range(c(rg.amax,rg.pot))

	plot(0, 0,
			 xlab = 'Return period (year)',
			 xlim = c(.5,6.5),
			 ylab = 'Flood quantiles (m^3/s)',
			 ylim = rg.qua, axes = FALSE)
	axis(2)
	axis(1, at = 1:6, c(2,5,10,20,50,100))


	if(out.amax){
	  arrows(1:6 - .1, qua.amax$lb, y1 = qua.amax$ub, code = 3, angle = 90, length = .05)
    points(1:6 - .1, qua.amax$qt, pch = 16, col = 'blue')
	}

  if(out.pot){
	  arrows(1:6 + .1, qua.pot$lb, y1 = qua.pot$ub, code = 3, angle = 90, length = .05)
	  points(1:6 + .1, qua.pot$qt, pch = 15, col = 'red')
  }

	legend('topleft', legend = c('AMAX','POT'),
				 col = c('blue','red'), pch = c(16,15))

	## Coefficient of variation ##
  if(out.amax){
	  cv.amax <- with(qua.amax, se/qt) * 100
	} else{
		cv.amax <- NULL
	}

	if(out.pot){
	  cv.pot <- with(qua.pot, se/qt) * 100
	} else{
		cv.pot <- NULL
	}

	cv.rg <- range(c(cv.pot, cv.amax))
	cv.rg <- cv.rg + c(-1,1) * 0.05 * diff(cv.rg)

	plot(0, 0,
			 xlab = 'Return period (year)',
			 xlim = c(.5,6.5),
			 ylab = 'Coefficient of variation (%)',
			 ylim = cv.rg, axes = FALSE)
	axis(2)
	axis(1, at = 1:6, c(2,5,10,20,50,100))

	if(out.amax){
  	lines(cv.amax)
	  points(cv.amax, col = 'blue', pch = 16)
		}

	if(out.pot){
  	lines(cv.pot)
	  points(cv.pot, col = 'red', pch = 15)
	}

	dev.off()

}

stopCluster(cl)

write(format(Sys.time() - t0), file = progress.file, append = TRUE)
write(format(Sys.time() - t0), file = error.file, append = TRUE)

##############################################
## Merge all the results obtained in parallel
##############################################

all.files <- list.files(DIR_CACHE, full.names = TRUE, pattern = '*.csv')
out <- lapply(all.files, read.csv)
out <- do.call(rbind,out)

outfile <- gzfile(file.path(DIR_OUT, 'atsite_flood_quantiles.csv.gz'))
## Uncomment to save new change
#write.csv(out, file = outfile, row.names = FALSE)


fail.files <- substr(list.files(DIR_CACHE, pattern = '*.png'), 1,7)
which(!(fail.files %in% GAUGEDSITES$station))


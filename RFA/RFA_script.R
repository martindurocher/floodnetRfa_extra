library(floodnetRfa)
library(CSHShydRology)
require(foreach)
require(doParallel)

rm(list = ls())

## Path of the HYDAT database
## Create a variable DB_HYDAT
source('config')

## Size of the bootstrap sample for L-moments methods
NSIM<- 1000

## Number of parallel processes
NCPU <- 4

## Resolution of the graphics
RES <- 360

## script and cache folder
DIR_OUT <- 'RFA'
DIR_CACHE <- 'cache/RFA'


################################################################################
## Monitoring file for parallel computing

if(!file.exists('cache'))
	dir.create('cache')

if(!file.exists(DIR_CACHE))
	dir.create(DIR_CACHE)

progress.file <- file.path(DIR_CACHE,"progress.txt")
error.file <- file.path(DIR_CACHE,"error.txt")

t0 <- Sys.time()
write(format(t0), file = progress.file)

## Set-up parallel computing
cl <- makeCluster(NCPU)
registerDoParallel(cl)

#############################################################
## Create a list of stationary station for each super regions
#############################################################

stations <- GAUGEDSITES$stations
sdistance <- SeasonDistanceData(DB_HYDAT, GAUGEDSITES$station)

## filter nonstationary sites
id.amax <- with(GAUGEDSITES, trend_mk >= 0.05 & trend_pt >= 0.05)
info.amax <- GAUGEDSITES[id.amax,]
supreg.amax <- with(info.amax, split(station, supreg_km12))

id.pot <- with(GAUGEDSITES, trend_lg >= 0.05 & trend_mx >= 0.05)
info.pot <- GAUGEDSITES[id.pot,]
supreg.pot <- with(info.pot, split(station, supreg_km12))

#############################################################
## Perform all RFA
#############################################################

statu <- foreach(ii = 1:4, #nrow(GAUGEDSITES),
								 .packages = c('floodnetRfa', 'CSHShydRology'),
								 .errorhandling = 'pass') %dopar%{


	site <- as.character(GAUGEDSITES$station[ii])
	supreg.id <- GAUGEDSITES$supreg_km12[ii]
	
	## allocate memory
	vars.ii <- data.frame(hamax = NA, namax = NA, erramax = TRUE,
											  hpot = NA, npot = NA, errpot = TRUE)

	write(paste0(ii, ' - Analyzing ', site),
				file = progress.file, append = TRUE)


	#######################
	## RFA - AMAX
	#######################

	## Extract seasonal distance among the station in the super region
	sreg.amax <- as.character(supreg.amax[[supreg.id]])
	sreg.amax <- sort(unique(c(site, sreg.amax)))
	sdist.amax <- sdistance[site, sreg.amax]

	## Extract annual maximums
	an <- AmaxData(DB_HYDAT, sreg.amax, distance = sdist.amax)

	## RFA
  amax <- try(FloodnetPool(an, target = site, nsim = NSIM,
  												 out.model = TRUE, verbose = FALSE),
  						silent = TRUE)
  
  
  ## Save results in DIR_CACHE
  if(class(amax) != 'try-error'){
		amax.file <- file.path(DIR_CACHE, paste0('rfAmax',site,'.csv'))
		write.csv(amax$qua, file = amax.file, row.names = FALSE)
		
		vars.ii$erramax <- FALSE
		vars.ii$hamax <- amax$fit$stat[1]
		vars.ii$namax <- nrow(amax$fit$lmom)

	} else {
		write(paste0('Error analyzing RFA-AMAX for station: ', site,']'),
					file = error.file, append = TRUE)
	}
  
 
  #######################
	## RFA - POT
	######################

  ## Extract seasonal distance 
  sreg.pot <- as.character(supreg.pot[[supreg.id]])
  sreg.pot <- sort(unique(c(site,sreg.pot)))
	sdist.pot <- sdistance[site,sreg.pot]
	
	## Extract stations info and peaks
	info <- GAUGEDSITES[GAUGEDSITES$station %in% sreg.pot, 
											c('station','auto','area')]

  xd <- DailyPeaksData(DB_HYDAT, info, distance = sdist.pot)

  ## RFA
  pot <- try(FloodnetPool(xd, target = site, nsim = NSIM, out.model = TRUE, 
  												verbose = FALSE), 
  					 silent = TRUE)


  ## Save results in DIR_CACHE
  if(class(pot) != 'try-error'){
		pot.file <- file.path(DIR_CACHE, paste0('rfPot',site,'.csv'))
		write.csv(pot$qua, file = pot.file, row.names = FALSE)

		vars.ii$errpot <- FALSE
		vars.ii$hpot <- pot$fit$stat[1]
		vars.ii$npot <- nrow(pot$fit$lmom)

	} else {
		write(paste0('Error analyzing RFA-POT for station: ', site,']'),
					file = error.file, append = TRUE)
	}
  
  ##############################
  ## Save heterogeneity measure statistics
  ##############################

  hetero.file <- file.path(DIR_CACHE,paste0('hetero', site, '.csv'))
  write.csv(vars.ii, file = hetero.file, row.names = FALSE)
  
  #######################
	## Summary Figure
	#######################

  Fd <- function(z){
  	ans <- z
  	colnames(ans) <- c('qt','se','lb','ub')
  	ans$cv <- ans$se/ans$qt *100

  	return(as.data.frame(ans))
  }

  ## Extract fitting statistics
  rgl <- list()
  if(!vars.ii$erramax){
    xa <- Fd(amax$qua)
    rgl$xa <- c(xa$lb,xa$ub)
  }
    
  if(!vars.ii$errpot){
    xp <- Fd(pot$qua)
    rgl$xp <- c(xp$lb,xp$ub)
  }
  
  
  fig.file <- paste0(site, '.png')
	png(file = file.path(DIR_CACHE, fig.file), height = 2*RES, width = 2*RES)

	par(mfrow = c(2,2))

	#### Pooling groups ####
	
	pid <- GAUGEDSITES$station %in% xd$site
	tid <- GAUGEDSITES$station %in% site
	scrd <- GAUGEDSITES[,c('season_x', 'season_y')]
	desc <- log(GAUGEDSITES[,c('area', 'map')])

	JulianPlot()
	title(main = 'Seasonal space')
	points(scrd, pch = '.')
	points(scrd[pid,], pch = 16, col = 'blue')
	points(scrd[tid,], pch = 10, col = 'red', cex = 2)

	plot(desc, pch = '.',
			 xlab = 'Drainage Area (log)', 
			 ylab = 'Mean annual precipitation (log)')
	title(main = 'Descriptor space')
	points(desc[pid,], pch = 16, col = 'blue')
	points(desc[tid,], pch = 10, col = 'red', cex = 2)

  
	#### Return level plot ####
	rg.qua <- range(unlist(rgl))
	rg.qua <- rg.qua + .02 * c(-1,1) * diff(rg.qua)

	prov <- GAUGEDSITES$prov[ii]
	
	plot(0, 0,
			 main = paste0(prov,' - ',site),
			 xlab = 'Return period (year)',
			 xlim = c(.5,6.5),
			 ylab = 'Flood quantiles (m^3/s)',
			 ylim = rg.qua, axes = FALSE)
	axis(2)
	axis(1, at = 1:6, c(2,5,10,20,50,100))

	Fcv <- function(z, pc, cl,  d){
		arrows(1:6 - d, z$lb, y1 = z$ub, code = 3, angle = 90, length = .05)
    points(1:6 - d, z$qt, pch = pc, col = cl)
	}
	
	if(!vars.ii$erramax) Fcv(xa,  16, 'blue',   -.1)
	if(!vars.ii$errpot)  Fcv(xp,  15, 'red',    .1)
	
	legend('topleft', legend = c('AMAX','POT'),
				 col = c('blue','red'), 
				 pch = c(16,15,17,18))

	#### Coefficient of variation ####
	
  if(!vars.ii$erramax){
	  cv.amax <- xa$cv
	} else {
		cv.amax <- NULL
	}

	if(!vars.ii$errpot){
	  cv.pot <- xp$cv
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

	if(!vars.ii$erramax){
  	lines(cv.amax)
	  points(cv.amax, col = 'blue', pch = 16)
		}

	if(!vars.ii$errpot){
  	lines(cv.pot)
	  points(cv.pot, col = 'red', pch = 15)
	}

	dev.off()

	return('exit normally')
}

stopCluster(cl)

write(format(Sys.time() - t0), file = progress.file, append = TRUE)

##############################################
## Merge all the results obtained in parallel
##############################################

all.files <- list.files(DIR_CACHE, full.names = TRUE, pattern = '*.csv')
rfa.files <- all.files[grep('rfAmax|rfPot',all.files)]
hetero.files <- all.files[grep('hetero',all.files)]

## rfa
out <- lapply(rfa.files, read.csv)
out <- do.call(rbind,out)
outfile <- gzfile(file.path(DIR_OUT, 'rfa_flood_quantiles.csv.gz'))

## Uncomment to save new change
#write.csv(out, file = outfile, row.names = FALSE)

## Hetero
out <- lapply(hetero.files, read.csv)
out <- do.call(rbind,out)[,c(1:2,4:5)]
outfile <- gzfile(file.path(DIR_OUT, 'heteo.csv.gz'))

## Uncomment to save new change
#write.csv(out, file = outfile, row.names = FALSE)

## verify for stations that have not been analyzed
cache.files <- substr(list.files(DIR_CACHE, pattern = '*.png'), 1,7)
which(!(GAUGEDSITES$station %in% cache.files))



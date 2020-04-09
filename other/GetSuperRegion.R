###############################################################################
#' Super regions
#'
#' Return the super region associated with a target site as suggested by the
#' dataset \code{gauged_targets.csv}.
#'
#' @param xd Imported data, i.e. \code{gauged_sites.csv}. 
#'  
#' @param target Station ID.
#'
#' @param type Type of frequency analysis: Annual maximum (\code{'amax'}) or
#'   peaks over threshold (\code{'pot'}).
#'
#' @param cluster Clustering thechniques. It corresponds to one column of
#'   \code{gauged_sites.csv}. Must be one of \code{'km6'}, \code{'km12'},
#'   \code{'hc6'} or \code{'hc12'}.
#'
#' @param trend.tol Significance level of the trend tests. This used to remove
#'   stations that may not have constant flood risk over time.
#'
#' @param thresh Threshold to return. It corresponds to a column of
#'   \code{gauged_sites.csv}. For example, \code{'auto'} or \code{'ppy200'}.
#'
#' @examples
#'
#' xd <- read.csv('gauged_sites.csv') 
#'  
#' GetSuperRegion('01AF009')
#
#' GetSuperRegion('01AD002', type = 'pot', thresh = 'ppy200')
#'
GetSuperRegion <-
	function(xd, target, 
			 type = 'amax',
			 cluster = 'km12',
			 trend.tol = 0.05,
			 thresh = 'auto'){

	## Find targets in the same super region
	sreg <- xd[,paste0('supreg_',cluster)]
	target.supreg <- sreg[xd$station == target]
	cond.supreg <- sreg == target.supreg


	if(type == 'amax'){
		## Identify the stationary targets
		cond.trend <- xd$trend_mk >= trend.tol &
									xd$trend_pt >= trend.tol

		## Return the intersection
		ans <- as.character(xd[cond.supreg & cond.trend, 'station'])

	} else if (type == 'pot'){
		## Identify the stationary targets
		cond.trend <- xd$trend_mx >= trend.tol &
									xd$trend_lg >= trend.tol

		## Return the intersection
		ans <- xd[cond.supreg & cond.trend, c('station', thresh, 'area')]
	}

	return(ans)
}

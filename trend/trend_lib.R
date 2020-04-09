#########################################################
#' Mann-Kendall's test with automatic selection of the lag
#' 
#' Function that perform the Mann-Kendall's test, and
#' evaluate the p-values using block bootstrap when there is 
#' autocorrelation.
#' The lag parameter is automatically identified when necessary.
#'
#' @param x Vector of data
#' @param lag.max Largest lag to verify
#' 
#' @return 
#' \itemize{
#'    \item tau : Statistic of the test.
#'    \item pvalue : P-value of the test.
#' }
#########################################################
AutoMK <- function(x, lag.max = 5){

	x <- as.numeric(x)
	nx <- length(x)

  ## Compute the PACF and
  ## Find the number of significant lag
  pa.ci <- 1.96 / sqrt(nx)
  pa <- pacf(x, na.action = na.pass, plot = FALSE)$acf[1:lag.max]
  pa.lag <- sum(abs(pa) > pa.ci)

  if(pa.lag == 0){
    mk <- Kendall::MannKendall(x)
    tau <- as.numeric(mk$tau)
    pval <- as.numeric(mk$sl)


  } else {

    ## evaluate p-value using block bootstrap
    MKtau<-function(z) Kendall::MannKendall(z)$tau
    mk.boot <- boot::tsboot(an$value, MKtau, R=2000, l=pa.lag+1, sim="fixed")
    tau <- mk.boot$t0
    pval <- mean(abs(mk.boot$t) > mk.boot$t0)
  }

  return(c(tau = tau, pvalue = pval))

}

#########################################################
#' Trend in exceedance rate using logistic regression
#' 
#' Verify that there is no trend in the number of peaks.
#' Logistic regression model with linear trend is applied and a t-test
#' the p-value of the t-test for the slope is returned.
#'
#' @param x Vector of covariate. Normally a time variable.
#' 
#' @param peaks Indices of the peaks in x
#' 
#' @return Minimal p-value
#########################################################

TrendLogis <- function(x, peaks){

  ## Trend in logistic regression
  xd <- data.frame(y = seq_along(x) %in% peaks, x = x)

  fit <- glm(y~x, xd, family = quasibinomial())
  ans <- summary(fit)$coefficient[2,4]

  ## Return minimal p-value
  return(ans)
}

## Verify using a F-test that trend model in the exceedance probability
## is better than the constant model
# TrendLogis <- function(x, peaks){
#
#   ## Trend in logistic regression
#   xbin <- data.frame(y = seq_along(x) %in% peaks, x = x)
#  
#   ## Constant model
#   fit0 <- glm(y~1, xbin, family = quasibinomial())
#
#   ## alternative models
#   fit <- vector('list', 3)
#   fit[[1]] <- glm(y~x, xbin, family = quasibinomial())
#   fit[[2]] <- glm(y~poly(x,2), xbin, family = quasibinomial())
#   fit[[3]] <- glm(y~poly(x,3), xbin, family = quasibinomial())
#
#   ## Evaluate p-values of the F-test
#   fun <- function(z) anova(fit0, z, test = 'F')[2,6]
#
#   ## Return minimal p-value
#   return(min(vapply(fit,fun, numeric(1)), na.rm = TRUE))
# }
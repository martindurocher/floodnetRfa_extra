###############################################################################
## objective:
## The function CSHShydRology::LmomDiag use polynomial approximation of the
## theoritical relation between the L-kurtosis and the L-skewness.
## The code below present the methodology to obtain the polynomial coefficients
###############################################################################

## Create a function that compute the theoritical kurtosis
library(lmomco)
Fkur <- function(s, d = 'nor'){
	k <- try(par2lmom(lmom2par(vec2lmom(c(1,1,s)), d))$ratio[4], silent = TRUE)

	if(is(k, 'try-error'))
		return(NA)
	else
		return(k)

}

## Evaluate the kurtosis for several skewness and distribution
tt <- seq(-.95,.95, 0.05)
dlst <- c('glo','gev','gno','pe3','gpa')

kur <- data.frame(matrix(0,length(tt),length(dlst)))
colnames(kur) <- dlst

for(d in dlst)
	kur[[d]] <- mapply(Fkur, tt, d)

## Add the limit points (1,1) and (-1,1)
kur <- data.frame(skew = tt, kur)
kur <- rbind(c(-1, rep(1,length(dlst))),
						 kur,
						 c( 1, rep(1,length(dlst))))

##############################################
## Vizualizing the results of the relation
##############################################

plot(gpa~skew, kur, type = 'l')
for(ii in seq(2,ncol(kur)-1))
	lines(kur$skew, kur[,ii], col = ii)

########################################################################
## Manually found the order of the polynomial that well fit the points
## The Pearson type III seems to be the most difficult to fit and appears
## relatively well with an order 10 polynomial
########################################################################
plot(glo~skew, kur)
lines(kur$skew, fitted(lm(glo ~ poly(skew, 10, raw = TRUE), kur)))

#################################################################
## Evaluate the magnitude of the errors
## The worst case scenario is 3e-4, which is satisfactory
#################################################################

Fres <- function(z) {
	f <- as.formula(paste(z,  "~ poly(skew, 10, raw = TRUE)"))
	res <- residuals(lm(f, kur))
	summary(abs(res))
}

for(s in c('gev','glo','gno','pe3', 'gpa')){
  print(s)
  print(Fres(s))
}
  

########################################################################
## Compute the polynomial coefficients to save
########################################################################

Fcoef <- function(z)	coef(lm(kur[,z] ~ poly(kur$skew, 10, raw = TRUE)))
out <- mapply(Fcoef, dlst)

## The previous output will be used to plot the l-moment ratio diagram

bmat <- read.table(text = "
             glo          gev           gno           pe3           gpa
0   1.666667e-01  0.107158699  1.225038e-01  1.226259e-01 -4.316343e-11
1   5.500064e-16  0.112615852  2.100695e-17  4.849455e-16  2.000000e-01
2   8.333333e-01  0.845179976  7.923906e-01  2.887205e-01  9.600000e-01
3  -1.085827e-14 -0.073172867 -7.453141e-15 -4.074033e-15 -1.919999e-01
4   2.658698e-14  0.005277763 -2.114303e-02  1.064608e+00  3.839996e-02
5   4.643463e-14 -0.066938343  4.843957e-14  1.743865e-14 -7.680604e-03
6  -9.309061e-14  0.100281501  2.770630e-01 -8.906645e-01  1.536177e-03
7  -7.209454e-14  0.087677507 -9.447091e-14 -3.055534e-14 -3.057862e-04
8   1.274086e-13 -0.132783946 -3.775294e-01  5.728072e-01  6.109048e-05
9   3.641720e-14 -0.060103832  5.399747e-14  1.810627e-14 -1.370841e-05
10 -5.915841e-14  0.074848480  2.066031e-01 -1.581181e-01  2.769001e-06
")

#################################################################
## Visualization of the fitting
#################################################################
sk <- seq(-1,1, len = 50)
xmat <- model.matrix(~poly(sk, 10, raw = TRUE))
pt <- xmat %*% as.matrix(bmat)

plot(sk, pt[,5], type = 'l')
for(ii in 1:5)
	lines(sk, pt[,ii],, col = ii+1)

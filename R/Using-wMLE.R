####################################################################
# This short script demonstrate the use of the weighted-MLE estimation
# technique within R (tested on R x64, 4.1.0).
#
# Always cite:
# Cousineau, D. (2009). Nearly unbiased estimators
#	for the three-parameter Weibull distribution
#	with greater efficiency than the iterative likelihood 
#	method. British Journal of Mathematical and 
#	Statistical Psychology, 62(1), 167-191.
####################################################################




# First: move to the folder containing the functions
# and source the function
source("C:/Users/DenisCousineau/Desktop/wMLE/R/wbl3.R")

# Second: initialize some data
suzuki = c(3.84,1.00,4.14,4.81,5.72,7.23,8.08,4.16,4.17,4.00,
		   4.42,3.58,3.92,4.73,5.42,5.09,5.59,3.67,5.76,6.34,
		   6.07,6.75,4.07,7.34,6.00,8.26,8.01,8.67,4.24,5.73,
		   5.50);

wbl3LogLikelihood(suzuki, c(2,1,-5))

# Third: run a regular MLE fit
wbl3MLE(suzuki)

# OR: run a weighted MLE fit
wbl3wMLE(suzuki)

# It is possible to specify some starting values such as
wbl3wMLE(suzuki, c(0.25, Inf, min(suzuki)-sd(suzuki)/2) )

#############################################################################
# Quick and dirty simulations; each takes roughly 10 minutes for 1000
# Note that in R, only a single starting value can be provided per parameter.
# In Mathematica, a bracket of starting values can be provided per parameter.
#############################################################################

# CASE A: wMLE
res1 <- c();
for( i in 1:5000 ) {
	X   <- rweibull(30, shape = 2, scale = 100)+300;
	sol <- wbl3wMLE(X, c(0.25, Inf, min(X)-sd(X)/2) )
	res1[i] <- sol$bestfitparams[1]
}
mean(res1)
# with theta =2,1,0    : 1.99
# with theta =2,100,300: 1.99;  2.033 was reported in 2009 using different heuristics

# CASE B: 2-parameter MLE with smallest datum removed based on MASS
library(MASS)
res2 <- c();
for( i in 1:5000 ) {
	X   <- sort(rweibull(30, shape = 2, scale = 100)+300);
	sol <- fitdistr((X-X[1])[-1], densfun="weibull");
	res2[i] <- sol$estimate[1]
}
mean(res2)
# with theta =2,1,0:     1.74; wrong of course because smallest data removed
# with theta =2,100,300: 1.75; same

# CASE C: regular MLE provided herein
res3 <- c();
for( i in 1:5000 ) {
	X   <- rweibull(30, shape = 2, scale = 100)+300;
	sol <- wbl3MLE(X, c(0.25, sd(X)/2, min(X)-sd(X)/2) )
	res3[i] <- sol$bestfitparams[1]
}
mean(res3)
# with theta =2,1,0:     1.96; good, owing to good heuristics
# with theta =2,100,300: 1.94; 1.829 was reported in 2009;





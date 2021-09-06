############################################################################
# Denis Cousineau, 2020, license CC 3.0
# version 0.1 09/05/2020.
# source this file to load the weights and have access to two functions
# 	wbl3MLE and wbl3wMLE
############################################################################


############################################################################
# subsidiary function: the log likelihood of the three-parameter Weibull distribution
wbl3LogLikelihood <- function(data, params) {
	sum( dweibull(data-params[3], shape= params[1], scale = params[2], log=TRUE ))
}
############################################################################




############################################################################
# The maximum likelihood fit of the 3-parameter Weibull distribution
# Input data: a vector of data
#       params, a three-item suggested starting values for the parameters:
#              gamma: the shape parameter, >0
#              beta: the scale parameter, >0
#              alpha: the shift (lower bound) parameter, alpha < min(data) 
#       options, that will be passed to the constrOptim function as is
#       
# Output: a list with the maximum loglikelihood fit of the data and the best-fitting parameters

wbl3MLE <- function(data, params = NULL, ...) {

    if (!(is.vector(data) & all(is.numeric(data)))) 
       stop('First argument must be a vector of numeric...');

    if  (!is.null(params)) {   		# explicit starting parameter values set by user
		if (length(params) != 3) 
			stop('params must contain precisely three parameters: c(gamma, beta, alpha)... Exiting.'); 
		gamma = params[1];
		beta  = params[2];
		alpha = params[3];
    } else {
		gamma = 2;              	# heuristic value
		beta  = sd(data)/2;     	# heuristic value
		alpha = 0.9 * min(data); 	# heuristic value
    }

    theta = c(gamma, beta, alpha); 	# put initial parameter values in an array
	objective <- function(params) -wbl3LogLikelihood(data, params)

    res <- constrOptim(
		theta,
		f  = objective,
        ui = diag(c(1, 1, -1)), ci = c(0, 0, -min(data)),
		method = "Nelder-Mead",
		...
    );

    # return likelihood to its natural direction below zero
    fit = -res$value;

    list(fit = fit, bestfitparams = res$par)
}
############################################################################














############################################################################
# The weighted maximum likelihood fit of the 3-parameter Weibull distribution
#
# The weights must be in a folder one level up called "weights".
# Limitations: the weights are currently pre-computed only up to n = 100.
#
# Based on Cousineau, D. (2009) Nearly unbiased estimators for the
#   three-parameter Weibull distribution with greater efficiency than the iterative 
#   likelihood method. British Journal of Mathematical and Statistical Psychology, 
#   62, 167–191. doi: 10.1348/000711007X270843
#
# Input data: a vector of data
#       params, a three-item suggested starting values for the parameters:
#              gamma: the shape parameter, >0
#              beta: the scale parameter, >0
#              alpha: the shift (lower bound) parameter, alpha < min(data)
#       options, that will be passed to the constrOptim function as is.
#       
# Output: a list with the w- MLE best-fitting parameters and 
#			the log-likelihood score (not a valid index of fit for wMLE)

wbl3wMLE <- function(data, params = NULL, ...) {

    ########################################################################
    ##  PART 0- Error checking
    ########################################################################

    if (!(is.vector(data) & all(is.numeric(data)))) 
       stop('First argument must be a vector of numeric...');
	n <- length(data);

    if (length(data) > dim(AllJ1)[1])
       stop(paste('The current implementation is limited to sample sizes from ',
			min(AllJ1[,1]), " to ", max(AllJ1[,1]), '... Exiting.', sep="") );

    if  (!is.null(params)) {   # explicit starting parameter values set by user
        if (length(params) != 3) 
            stop('params must contain precisely three parameters: c(gamma, beta, alpha) (even if beta is unused)... Exiting.'); 
        gamma = params[1];
        # beta  = params[2];     # unused in w-MLE
        alpha = params[3];
    } else {
        gamma = 2;               # heuristic value: 1 is exponential; 10 is degenerated
        # beta  = 2 * sd(data);  # unused in w-MLE
        alpha = 0.9 * min(data); # heuristic value
    }

    # cat("End of PART 0\n");
    # cat("   dimensions: ", n, "\n");
    # cat("   heuristics: ", gamma, ", ", alpha, "\n");


    ########################################################################
    ##  PART 1- loading embedded functions here...
    ########################################################################
    # these functions are defined in Cousineau (2009)

    fctGamma <- function(x, n, gamma, alpha, W2) {
        W2/gamma + sum( log(x - alpha ) )/n -
               sum((log(x - alpha ) * (x -alpha)^gamma ) / sum((x -alpha )^gamma )) 
    }

    fctBeta <- function(x, n, gamma, alpha, W1) {
         (sum((x - alpha )^gamma )/(n * W1))^(1/gamma)
    }

    fctAlpha <- function(x, n, gamma, alpha, W3) {
         -W3 + ( sum(1 / (x - alpha )) * sum( (x - alpha)^gamma)  / sum( (x - alpha)^(gamma-1)) ) / n
    }

    # cat("End of PART 1\n");
    # cat("   test fctGamma:  ", fctGamma(c(1,2,3),3,2,0,5), "\n");
    # cat("   test fctBeta:   ", fctBeta (c(1,2,3),3,2,0,5), "\n");
    # cat("   test ftcAlpha:  ", fctAlpha(c(1,2,3),3,2,0,5), "\n");


    ########################################################################
    ##  PART 2- loading precomputed weights...
    ########################################################################

    # these two depends only on sample size n and thus are constants herein
    J1 <- AllJ1[AllJ1[1] == n, 2];
    J2 <- AllJ2[AllJ2[1] == n, 2];

    # the third one depends on gamma as well and so must be dynamic
    sub <- AllJ3[AllJ3[1] == n, , ];
    # this function does linear interpolation
	J3 <- function(gamma) {
		if (gamma<0.1) {gamma = 0.1};
		if (gamma>5.0) {gamma = 5.0};
		pos = floor(gamma / wMLEresolution )+1; # gamma of 0 are in position 1
		low = sub[pos-1, ];
		hig = sub[pos,   ];
		# returns interpolation
		as.numeric(((gamma - low[2])/wMLEresolution * (hig[3]-low[3]) + low[3]))
	}

    # cat("End of PART 2\n");
    # cat("   Weights: ", J1, J2, J3(2.66131), "\n" );


    ########################################################################
    ##  PART 3- Set up the search parameters
    ########################################################################

	# nothing to do with R

    theta <- c(gamma, alpha);                # put initial parameter values in an array

    # cat("End of PART 3\n");
    # cat("   theta: ", theta, "\n" );


    ########################################################################
    ##  PART 4- Run the search and conclude
    ########################################################################

    objective <- function(theta) {
        (fctGamma(data, n, theta[1], theta[2], J2))^2 + (fctAlpha(data, n, theta[1], theta[2], J3(theta[1]) ))^2
    }

    # run the search over the two parameters gamma and alpha
    res <- constrOptim(
		theta,
		f      = objective,
        ui     = diag(c(1, -1)), ci = c(0, -min(data)),
		method = "Nelder-Mead",
		...
    );
    # cat("value and steps ", res$value, res$counts[1], "\n");

    # unpack the two parameters and solve for the third
    gamma <- res$par[1];
    alpha <- res$par[2];
    beta  <- fctBeta(data, n, gamma, alpha, J1);

    # this is it
    # cat("End of PART 4\n");
    # cat("   gamma:  ", gamma, "\n" );
    # cat("   beta:   ", beta, "\n"  );
    # cat("   alpha:  ", alpha, "\n" );
    # cat("   sol:    ", objective(res$par), "\n"); # always 0

    bestfitparams = c(gamma, beta, alpha);
    fit = wbl3LogLikelihood(data, bestfitparams); #not a real fit index

	list(fit = fit, bestfitparams = bestfitparams)
}


########################################################################
##  All done. What follows are the lookup tables
########################################################################

initwts <- function(path) {
    # a matrix of pairs {n, J1(n)} for n going from 1 to 100
    AllJ1 <<- read.table(paste(path,'../weigths/J1.tsv',sep="/"),sep="\t");

    # a matrix of pairs {n, J2(n)} for n going from 1 to 100
    AllJ2 <<- read.table(paste(path,'../weigths/J2.tsv',sep="/"),sep="\t");

    # a matrix of triplets {n, gamma, J3(n, gamma)} for n going from 1 to 100
    AllJ3 <<- read.table(paste(path,'../weigths/J3.tsv',sep="/"),sep="\t");

    # the size of the mesh in the lookuptable
    wMLEresolution <<- AllJ3[2,2]-AllJ3[1,2];

}

# launch initialization
initwts("C:/Users/DenisCousineau/Desktop/wMLE/R")

cat("All good to go. Weights are loaded for sample sizes up to ", dim(AllJ1)[1], "\n")

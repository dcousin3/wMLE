This folder contains the weights to be used with the 
w-MLE estimator.

The ones used in the scripts are J1, J2, J3, that is, the
bias correcting weights in the median sense. Are also 
included G1, G2 and G3 (in the geometric mean sense just 
because I had them under the hand...)

The weights must accompany the scripts, and be located
in a folder called "weights" located one level up the 
folder containing the script.

As of now, the weights are for sample sizes up to 100
only (and regarding the third file, for shape parameter
up to 5 only by step of 0.10). I may proceed to expand
(to larger n) and refine (with smaller steps) these 
weights in a near future.

The Mathematica program _GeneratingTheWeights.nb is used
to generate the weights in case more scenarios are needed. 

Always cite:
Cousineau, D. (2009). Nearly unbiased estimators
	for the three-parameter Weibull distribution
	with greater efficiency than the iterative likelihood 
	method. British Journal of Mathematical and 
	Statistical Psychology, 62(1), 167-191.

Let me know if this is usefull.

Denis
denis.cousineau@uottawa.ca


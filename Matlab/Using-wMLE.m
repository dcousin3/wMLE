% This short script demonstrate the use of the weighted-MLE estimation
% technique within Matlab (tested on R2018a).

% Always cite:
% Cousineau, D. (2009). Nearly unbiased estimators
%	for the three-parameter Weibull distribution
%	with greater efficiency than the iterative likelihood 
%	method. British Journal of Mathematical and 
%	Statistical Psychology, 62(1), 167-191.

% To use the fitting procedure with Matlab
% First: move to the folder containing the functions

% Second: initialize some data
suzuki = [3.84,1.00,4.14,4.81,5.72,7.23,8.08,4.16,4.17,4.00,4.42,3.58,3.92,4.73,5.42,5.09,5.59,3.67,5.76,6.34,6.07,6.75,4.07,7.34,6.00,8.26,8.01,8.67,4.24,5.73,5.50]

% Third: run a regular mle fit
[fit, params] = wbl3MLE(suzuki, [2.08, 3.09, 2.80])

% OR: run a weighted mle fit
[fit, params] = wbl3wMLE(suzuki)

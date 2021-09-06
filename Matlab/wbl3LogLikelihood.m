function logL = wbl3LogLikelihood(data, params)

% the loklikelihood of the 3-parameter Weibull distribution
% Denis Cousineau, 2020, license CC 3.0
% version 0.1 09/05/2020.

% Input data: a vector of data
%       params, a three-item vector containing the parameters:
%              gamma: the shape parameter, >0
%              beta: the scale parameter, >0
%              alpha: the shift (lower bound) parameter, alpha < min(data)
% Output logL: the loglikelihood of the data given the parameters

if (nargin ~= 2)
  error('There must be data followed by the parameter vector... Exiting.');
end
[n,m]=size(data);
if (n > 1) && (m > 1)
   error('First argument must be a datum or a vector of data... Exiting.');
end
if (n == 1) && (m > 1)       % case of a row vector
   data = data';
   n = m;
end
p = length(params);
if (not(isvector(params)) || (p ~= 3) )
  error('params must contain precisely three parameter gamma, beta, alpha... Exiting.'); 
end

ps   = wbl3PDF(data, params) + realmin;
logL = sum(log(ps)); 

return %logL

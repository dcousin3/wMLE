function p = wbl3PDF(data, params)

% the PDF of the 3-parameter Weibull distribution
% Denis Cousineau, 2020, license CC 3.0
% version 0.1 09/05/2020.

% Input data: a scalar or a vector of data
%       params: a vector of three parameters in that order:
%           gamma: the shape parameter, >0
%           beta: the scale parameter, >0
%           alpha: the shift (lower bound) parameter, alpha < min(data)
% Output p: a scalar or a vector (depending on data) of density

if (nargin ~= 2)
  error('There must be a datum followed by the parameters vector... Exiting.');
end

[n,m] = size(data);
if (n > 1) && (m > 1)
   error('First argument must be a datum or a vector of data... Exiting.');
end
if (n == 1) && (m > 1)     %case of a row vector of data
  data = data';
end
gamma = params(1);
beta  = params(2);
alpha = params(3);

if (gamma < 0) error ('The shape parameter must be positive... Exiting.'); end
if (beta < 0)  error ('The scale parameter must be positive... Exiting.'); end
if min(data) <= alpha     % case where alpha is smaller than the shift  
  error('At least one datum is below the lower bound parameter... Exiting.');
end

p = (beta^(-gamma) .* (exp(-((data - alpha)./beta).^gamma)) .* gamma .*(data - alpha).^(gamma-1)) ;

return %p

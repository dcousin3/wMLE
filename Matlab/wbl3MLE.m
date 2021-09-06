function [fit, bestfitparams] = wbl3MLE(data, params, options)

    % the maximum likelihood fit of the 3-parameter Weibull distribution
    % Denis Cousineau, 2020, license CC 3.0
    % version 0.1 09/05/2020.

    % Input data: a scalar or a vector of data
    %       params, a three-item suggested starting values for the parameters:
    %              gamma: the shape parameter, >0
    %              beta: the scale parameter, >0
    %              alpha: the shift (lower bound) parameter, alpha < min(data) 
    %       options, a three-item vector with tolerance on X, 
    %              tolerance on Fct,
    %              and maxIter
    %       
    % Output f: the maximum loglikelihood fit of the data and the best-fitting parameters

    [n,m]=size(data);
    if (n > 1) && (m > 1)
       error('First argument must be a datum or a vector of data...');
    end
    if  n == 1        	                   %case of a row vector of data
      data = data';
    end

    if  (nargin > 1) && ~isempty(params)   % explicit starting parameter values set by user
      p = size(params);
      if (p ~= 3) 
        error('params must contain precisely three parameters: gamma, beta, alpha... Exiting.'); 
      end
      gamma = params(1);
      beta  = params(2);
      alpha = params(3);
    else
      gamma = 2;               % heuristic value
      beta  = std(data)/2;     % heuristic value
      alpha = 0.9 * min(data); % heuristic value
    end

    % set the search parameters
    if  (nargin > 2) && ~isempty(options)    % user-defined options values passed to fmins
      opts = options(1:3);             % termination, function tolerances and maximum number of iterations
    else
      % default values for termination, function tolerances and default max number of iterations
      opts = [ 1.e-4, 1.e-4, 200 * length(data) ]; 
    end

    optionsfmin = optimset( "Display", "off", ...
      'TolX',opts(1),'TolFun',opts(2), ...
      'MaxFunEvals',opts(3),'MaxIter',opts(3) ...
    );

    pinit = [gamma beta alpha];                % put initial parameter values in an array

    warning ("off");
    [bestfitparams, fit, cvg, outp] = fmincon(@(params) -wbl3LogLikelihood(data, params), ...
                             pinit, [],[],[],[], ...
                             [0 0 -Inf]',[Inf, Inf, min(data)-0.00001]',[],optionsfmin ...
                            );
    warning ("on");

    % return likelihood to its natural direction below zero
    fit = -fit;

    return %[fit, bestfitparams]

end
function [fit, bestfitparams] = wbl3wMLE(data, params, options)

    % the weighted maximum likelihood fit of the 3-parameter Weibull distribution
    % Denis Cousineau, 2020, license CC 3.0
    % version 0.1 15/05/2020.
    % version 0.2 30/08/2021: uses a simplex instead of fmincon
    %
    % Based on Cousineau, D. (2009) Nearly unbiased estimators for the
    %   three-parameter Weibull distribution with greater efficiency than the iterative 
    %   likelihood method. British Journal of Mathematical and Statistical Psychology, 
    %   62, 167–191. doi: 10.1348/000711007X270843

    % Input data: a scalar or a vector of data
    %       params, a three-item suggested starting values for the parameters:
    %              gamma: the shape parameter, >0
    %              beta: the scale parameter, >0
    %              alpha: the shift (lower bound) parameter, alpha < min(data)
    %       options, a three-item vector with tolerance on X, 
    %              tolerance on Fct,
    %              and maxIter
    %       
    % Output f: the w-MLE best-fitting parameters
    % Limitations: the weights are currently pre-computed only up to n = 100.

    % This function is self-contained, i.e., all the required sub-functions
    % are embedded within it.
    % The weights must be in a folder one level up called "weights"


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%  PART 0- Error checking
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [n,m]=size(data);
    if (n > 1) && (m > 1)
       error('First argument must be a datum or a vector of data... Exiting.');
    end
    if  n == 1        	                   %case of a row vector of data
      data = data';
      n = m;
    end

    if (n > 64)
       error ('The current implementation is limited to n from 1 t 64... Exiting.');
    end

    if  (nargin > 1) && ~isempty(params)   % explicit starting parameter values set by user
        p = size(params);
        if (p ~= 3) 
            error('params must contain precisely three parameters: gamma, beta, alpha... Exiting.'); 
        end
        gamma = params(1);
        % beta  = params(2);     % unused in w-MLE
        alpha = params(3);
    else
        gamma = 2;               % heuristic value
        % beta  = 2 * std(data); % unused in w-MLE
        alpha = 0.9 * min(data); % heuristic value
    end

    % fprintf("End of PART 0\n");
    % fprintf("   dimensions: %d %d\n", n, m);
    % fprintf("   heuristics: %f %f %f\n", gamma, beta, alpha);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%  PART 1- loading embedded functions here...
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % these functions are defined in Cousineau (2009)

    function f = fctGamma(x, n, gamma, alpha, W2)
         f =  W2/gamma + sum( log(x - alpha ) )/n ...
              - sum((log(x - alpha ) .* (x -alpha).^gamma ) / sum((x -alpha ).^gamma )) ;
    end

    function f = fctBeta(x, n, gamma, alpha, W1)     
         f =  (sum((x - alpha ).^gamma )/(n * W1)).^(1/gamma);
    end

    function f = fctAlpha(x, n, gamma, alpha, W3)
         f =  -W3 + ( sum(1 ./ (x - alpha )) * sum( (x - alpha).^gamma)  / sum( (x - alpha).^(gamma-1)) ) / n;
    end    

    % fprintf("End of PART 1\n");
    % fprintf("   test fctGamma:  %f\n", fctGamma([1,2,3],3,2,0,5));
    % fprintf("   test fctBeta:   %f\n", fctBeta ([1,2,3],3,2,0,5));
    % fprintf("   test ftcAlpha:  %f\n", fctAlpha([1,2,3],3,2,0,5));


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%  PART 2- loading precomputed weights...
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    global AllJ1; 
    global AllJ2;
    global AllJ3;

    [resolution, ming, maxg] = initwts();

    % these two depends only on sample size n
     J1 = AllJ1(n, 2);
     J2 = AllJ2(n, 2);
 
     % the third one depends on gamma as well and so must be dynamic
     idx = any(AllJ3 == n, 2);
     sub = AllJ3(idx, 1:3);
     % this function does linear interpolation
     function f = J3(gamma) 
         if (gamma > maxg-resolution) gamma = maxg-2*resolution; end
         pos = floor(gamma / resolution )+1; % gamma of 0 are in position 1
         low = sub(pos, 1:3);
         hig = sub(pos+1,   1:3);
         f   = ((gamma - low(2))/resolution * (hig(3)-low(3)) + low(3));
     end

     fprintf("End of PART 2\n");
     fprintf("   resolution/maxg/ming %f / %f / %f\n", resolution, maxg, ming);
     fprintf("   Weights: %f %f %f\n", J1, J2, J3(2.66131   ) );


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%  PART 3- Set up the search parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % set the search parameters
    if  (nargin > 2) && ~isempty(options)    % user-defined options values passed to fmins
      opts = options(1:3);             % termination, function tolerances and maximum number of iterations
    else
      % default values for termination, function tolerances and default max number of iterations
      opts = [ 1.e-10, 1.e-10, 500 * length(data) ]; 
    end

    optionsfmin = optimset( "Display", "off", ...
      'TolX',opts(1),'TolFun',opts(2), ...
      'MaxFunEvals',opts(3),'MaxIter',opts(3) ...
    );

    pinit = [gamma alpha];                % put initial parameter values in an array

    % fprintf("End of PART 3\n");
    % fprintf("   pinit: %f %f \n", pinit );


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%  PART 4- Run the search and conclude
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function fit = objfct(data, p) 
        fit = (fctGamma(data, length(data), p(1), p(2), J2))^2 + (fctAlpha(data, n, p(1), p(2), J3(p(1)) ))^2;
    end

    % run the search for the first two parameters
    warning ("off");
    [bestfitparams, ~, ~, ~] = fmincon(@(p) objfct(data, p), ...
                                       pinit, [], [], [], [], ...
                                       [0.0 -Inf]', [5 min(data)]', [], optionsfmin ...
                               );
    warning ("on");
    % simplex with boundaries; step size alpha was reduced
%    dbg=0;
%    [bestfitparams, fval, ncall] = simplexmin( ...
%        @(p) objfct(data, p), ...
%        pinit, ...
%        optionsfmin, ...
%        [],[],[],[],[0.0 -Inf]', [5.0 min(data)]',dbg ...
%    );
    % fprintf("value and steps %f %d\n", fval, ncall);

    % unpack the two parameters and solve for the third
    gamma = bestfitparams(1);
    alpha = bestfitparams(2);
    beta  = fctBeta(data, n, gamma, alpha, J1);

    % this is it
    % fprintf("End of PART 4\n");
    % fprintf("   gamma:  %f\n", gamma );
    % fprintf("   beta:   %f\n", beta  );
    % fprintf("   alpha:  %f\n", alpha );
    % fprintf("   sol:    %f\n", objfct(data, bestfitparams)); % always 0

    bestfitparams = [gamma, beta, alpha];
    fit = wbl3LogLikelihood(data, bestfitparams); %not a real fit index

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  All done. What follows are the lookup tables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [resolution, ming, maxg] = initwts()

    global AllJ1; 
    global AllJ2;
    global AllJ3;

    % a matrix of pairs {n, J1(n)} for going from 1 to 64
    AllJ1 = cell2mat(textscan(fopen('..\weigths\J1.tsv'),'%f %f'));

    % a matrix of pairs {n, J2(n)} for going from 1 to 64
    AllJ2 = cell2mat(textscan(fopen('..\weigths\J2.tsv'),'%f %f'));

    % a matrix of triplets {n, gamma, J3(n, gamma)} for going from 1 to 64
    AllJ3 = cell2mat(textscan(fopen('..\weigths\J3.tsv'),'%f %f %f'));

    % the size of the mesh in the lookuptable
    resolution = AllJ3(2,2)-AllJ3(1,2);

    temp = max(AllJ3);
    maxg = temp(2);
    temp = min(AllJ3);
    ming = temp(2);

end

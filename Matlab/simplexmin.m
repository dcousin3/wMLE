function [x, fv, nfun] = simplexmin(fun,x,opt,Ane,Bne,Aeq,Beq,lb,ub,dbg)
% gradient minimizer using sequential simplex - straight simplex solver
% function [x, fv, nfun] = simplexmin(fun,x,opt,Ane,Bne,Aeq,Beq,lb,ub,dbg)
%
% inputs  10 - 8 optional
% fun     function to evaluate                          class char
% x       initial conditions                            class real
% opt     options (fminsearch)                          class struct
% Ane     inequality A matrix                           class real
% Bne     inequality B vector                           class real
% Aeq     equality A matrix                             class real
% Beq     equality B vector                             class real
% lb      lower bound on states                         class real
% ub      upper bound on states                         class real
% dbg     debug flag                                    class int
%
% outputs 3
% x       final states                                  class real
% fv      final cost                                    class real
% nfun    # of function evals                           class int
%
% michael arant - March 1995
% modified for explicit bound - July 2018
% added nonlinear inequality / equality - Oct 2018
%
% inequality:
% A*x <= B
% or
% Ane = {@(x)(x(1)^2+x(2)^2); @(x)(x(1)*x(2)^2)}
% Bne = {0; 3}
%
% equality:
% A*x = B
% or
% Aeq = {@(x)(x(1)^2+x(2)^2); @(x)(x(1)*x(2)^2)}
% Beq = {0; 3}
%
% inequality can be used in place of bounds if softer edges are desired:
% Ane = [-eye(numel(x)); eye(numel(x))];
% Bne = [-1 * lower_bounds; upper bounds];
%
% simplex approach based on 
% sequential simplex method - "Elements of Structural Optimization"
% R. Haftka / Z. Gurdal - pg 124 - 127
%
% examples
% [x, fv, nfun] = simplex(@(x)(x(1)^2+x(2)^2),[-3 4],[],[1 0; 0 1],[-1; -1],[],[],[],[],1)
% [x, fv, nfun] = simplex(@(x)(x(1)^2+x(2)^2),[-3 4],[],[1 0],[-1],[0 1],[.2],[],[],1)
% [x, fv, nfun] = simplex(@(x)(x(1)^2+x(2)^2),[-3 4],[],...
%	{@(x)(x(1)); @(x)(x(2))},[-1; -1],[],[],[],[],1)
% [x, fv, nfun] = simplex(@(x)(x(1)^2+x(2)^2),[-3 4],[],[],[], ...
%	{@(x)(x(1)); @(x)(x(2))},[-1; -1],[],[],1)
%% dbug
if ~exist('dbg','var') || isempty(dbg)
	dbg = 0;
end
%% inputs
% function test
if ~exist('fun','var') || isempty(fun)
	help(mfilename); x = []; fv = []; return
end
% ic test - if no ic given, abort
if ~exist('x','var') || isempty(x)
	help(mfilename); x = []; fv = []; return
end
% force column vector form
x = reshape(x,numel(x),1); len = numel(x);
% options - solver options if desired
if ~exist('opt','var') || isempty(opt) || ~isstruct(opt)
	% init the structure
	clear opt;
	opt = struct;
end
% set max itterations?
if ~isfield('MaxIter',opt)
	opt.MaxIter = 1e4 * numel(x);
end
% set function tolerance?
if ~isfield('TolFun',opt)
	opt.TolFun = 1e-4;
end
% set step tolerance?
if ~isfield('TolX',opt)
	opt.TolX = 1e-4;
end
% inequality test - if inequality not passed, set to null
if ~exist('Ane','var') || isempty(Ane)
	Ane = [];
end
if ~exist('Bne','var') || isempty(Bne)
	Bne = [];
end
if size(Ane,1) ~= size(Bne,1)
	error('Ane and Bne are different sizes')
end
% equality test - if equality not passed, set to null
if ~exist('Aeq','var') || isempty(Aeq)
	Aeq = [];
end
if ~exist('Beq','var') || isempty(Beq)
	Beq = [];
end
if size(Aeq,1) ~= size(Beq,1)
	error('Aeq and Beq are different sizes')
end
% bounds - if bounds not passed, set bound flag to 0
% upper bound test
if exist('ub','var') && ~isempty(ub)
	% force column vector
	ub = reshape(ub,len,1);
	% check for lower bound
	if ~exist('lb','var') || isempty(lb)
		lb = -Inf(size(x));
	end
end
% lower bound test
if ~exist('lb','var') || isempty(lb)
	% no bounds - set flag to zero
	bf = 0; ub = []; lb = [];
else
	% force column vector
	lb = reshape(lb,len,1); bf = 1;
	% crsoo check upper bound
	if ~exist('ub','var') || isempty(ub)
		ub = Inf(size(x));
	end
end
% set bound offset percentatge (this is the "step back" from the bound edge)
bo  = 0.01;
% check for bound crossing (lower bound above upper bound)
if exist('lb','var') && ~isempty(lb)
	if sum(lb > ub); error('lower bound > upper bound'); end
	if sum(lb > x); error('IC violates bounds'); end
	if sum(ub < x); error('IC violates bounds'); end
end
% test ic - push ic into bounds if initial IC was out of bound
if bf; x = bound(x,lb,ub,1,bo); end
%% function test
fun = fcnchk(fun);
%% default terms
% alpha is a simplex step
% beta is an expansion
% gamma is a contraction
% eta is a new simplex generation
% gs is the initial gradient step size
% alpha = 1; beta = alpha * 2; gamma = alpha / 2; eta = 0.5; gs = 0.05;
alpha = 0.2;        %%% important: reduce this step size 
beta = alpha * 2; beta=2;
gamma = alpha / 2; gamma =0.5;
eta = 0.5; 
gs = 0.5;
% Gao / Han
% nd = numel(x);
% alpha = 1; beta = 1 + 2/nd; gamma = 0.75 + 1/(2 * nd); eta = 1 - 1/nd; gs = 0.05;
% [alpha beta gamma eta]
%% start the solver
% initial function evaluation
fv = efun(fun,x,Ane,Bne,Aeq,Beq); nfun = 1;
if dbg
	fprintf('Initial function: %g	/	IC:	',fv); fprintf('%g,	',x);
	fprintf('\n');
end
% solve
[fv, x, nfun] = simplex_loop(fv,x,nfun,fun,opt, ...
			len,bo,bf,lb,ub,Ane,Bne,Aeq,Beq, ...
			alpha,beta,gamma,eta,gs,dbg);
%% simplex function
function [fv, x, nfun] = simplex_loop(fv,x,nfun,fun,opt,len,bo,bf,lb,ub,Ane, ...
	Bne,Aeq,Beq,alpha,beta,eta,gamma,gs,dbg)
% set initial function evaluation and point to infinity (force loop to start)
xp = Inf(size(x)); fp = Inf;
while max(abs(xp - x)) > opt.TolX && (fp - fv) > opt.TolFun ...
		&& nfun < opt.MaxIter
	
	% record last vector and function value
	xp = x; fp = fv;
	
	% size gradient matrix and function vector
	xg = x * ones(1,len+1); fg = NaN(1,len+1);
	
	% build gradients
	temp = len+1:len+1:(len+1)*len;
	xg(temp) = xg(temp) * (1 + gs);
	% zero test - X * 0 = 0....  so any 0 values are set to gs.
	xg(temp(~xg(temp))) = gs;
	
	% bound test if needed
	if bf; xg = bound(xg,lb,ub,len+1,bo); end
	
	% gradient evaluation - use a loop as the function may not support vectors
	for ii = 1:len+1
		fg(ii) = efun(fun,xg(:,ii),Ane,Bne,Aeq,Beq);
	end
	% number of fucntion evals
	nfun = nfun + len + 1;
	
	% mark initial max and min
	[~, minv] = min(fg);
	[~, maxv] = max(fg);
	% solve simplex - this is the clasical simplex problem
	ftol = Inf;
	while nfun < opt.MaxIter && ftol > opt.TolFun
		
		% mean point in design space NOT including the worst point
		xm = mean(xg(:,~((1:len+1) == maxv)),2);
		% bound test
		if bf; xm = bound(xm,lb,ub,1,bo); end
		% function eval
		fm = efun(fun,xm,Ane,Bne,Aeq,Beq); nfun = nfun + 1;
		
		% reflection
		xr = xm + alpha*(xm - xg(:,maxv));
		% bound test
		if bf; xr = bound(xr,lb,ub,1,bo); end
		% function eval
		fr = efun(fun,xr,Ane,Bne,Aeq,Beq); nfun = nfun + 1;
		
		% test step - right direction? reflection is better than
		% previous best
		if fr < fg(minv)
			% good step - expansion
			xe = xm + beta*(xr - xm);
			% bound test
			if bf; xe = bound(xe,lb,ub,1,bo); end
			% function eval
			fe = efun(fun,xe,Ane,Bne,Aeq,Beq); nfun = nfun + 1;
			
			% step still in the right difrection (expansion improve result)?
			if fe < fr
				% expansion step
				% update and repeat - replace max value with reflection
				xg(:,maxv) = xe; fg(maxv) = fe;
			else
				% reflection step
				% replace worst point (reflection was good, expansion was not)
				xg(:,maxv) = xr; fg(maxv) = fr;
			end
		else
			% reflection was not better than previous best guess
			% reflectin greater than second worst value?
			[~, sq] = sort(fg);
			if fr < fg(sq(end-1))
				% value was better than the second worst value.  repalce / repeat
				xg(:,maxv) = xr; fg(maxv) = fr;
			else
				% reflection better than the last worst value?
				if fr < fg(maxv)
					% value was better than the worst value. 
					% contraction - small step to re-test current space
					xc = xm - gamma*(xg(:,maxv) - xm);
				else
					% reverse direction 
					% contraction - small step to re-test current space
					xc = xm + gamma*(xg(:,maxv) - xm);
				end
				% bound test
				if bf; xc = bound(xc,lb,ub,1,bo); end
				% function eval
				fc = efun(fun,xc,Ane,Bne,Aeq,Beq); nfun = nfun + 1;
				
				% contraction was a better result?
				if fc > fg(:,maxv)
					% not a better result - build new simplex and retry 
					% shrinkage - new nominal point based on current best point
					xg = xg + eta*(xg(:,minv)*ones(1,len+1) - xg);
					% bound check
					if bf; xg = bound(xg,lb,ub,len+1,bo); end
					
					% new gradient fucntion evaluation
					for ii = 1:len+1
						fg(ii) = efun(fun,xg(:,ii),Ane,Bne,Aeq,Beq);
					end
					% update function count
					nfun = nfun + len + 1;
				else
					% contraction was a better result.  replace worst point
					xg(:,maxv) = xc; fg(maxv) = fc;
				end
			end
		end
		
		% max / min location for new simplex
		[~, maxv] = max(fg);
		[~, minv] = min(fg);
		
		% calcualte function tolerance for exiting - eqn 4.2.8
		ftol = sqrt(sum((fg - fm).^2) / (len+1));
		
		% dbg?
		if dbg > 1
			fprintf('min val / num fun =	%g	/	%g	',fg(minv),nfun);
			fprintf('Best point	');
			fprintf('%g,	',xg(:,minv)); fprintf('\n');
		end
	end
	
	% extract the resutls
	x = xg(:,minv); fv = fg(minv);
	
	if dbg
		fprintf('min val / num fun = %g	/	%g	',fv,nfun); 
		fprintf('Best point	'); fprintf('%g,	',x);
		fprintf('\n');
	end
end
%% function evaluation
function [f] = efun(fun,x,Ane,Bne,Aeq,Beq)
% disp values
% fprintf('%g	',x);
% evaluate function
f = fun(x);
% inequality test
if ~isempty(Ane)
	% assess the inequality - if r > 0 a violation exists
	if isa(Ane,'double')
		r = Ane*x - Bne; r(isnan(r)) = 0; r(r <= 0) = 0;
	else
		% evaluate nonlinear inequality
		r = NaN * Bne;
		for ii = numel(r):-1:1
			r(ii) = Ane{ii}(x) - Bne(ii);
		end
		r(isnan(r)) = 0; r(r <= 0) = 0;
	end
	% offset for better convergence
% 	fprintf('%1.5f	x = ',f); fprintf('%1.5f	',x); fprintf(' / ');
% 	fprintf('%1.5f	',r); fprintf(' / ');
	r = abs(r) + 1;
% 	disp([x' r'])
	% cost function
	f = f + (max(abs(f),1))*(10^(sum(r.^2)) - 10^numel(r));
% 	fprintf('%1.f\n',f);
end
% equality test
if ~isempty(Aeq)
	% asses the equality - if r ~= 0 a violation exits
	if isa(Aeq,'double')
		r = Aeq*x - Beq; r(isnan(r)) = 0;
	else
		% evaluate nonlinear equality
		r = NaN * Beq;
		for ii = 1:numel(r)
			r(ii) = Aeq{ii}(x) - Beq(ii);
		end
		r(isnan(r)) = 0;
	end
% 	fprintf('%1.5f	',f); fprintf('%1.5f	',x); fprintf(' / ');
% 	fprintf('%1.5f	',r); fprintf(' / ');
	% offset for better convergence
	r = abs(r) + 1;
% 	disp([x' r'])
	% cost function
	f = f + (max(abs(f),1))*(10^(sum(r.^2)) - 10^numel(r));
% 	fprintf('%g\n',f);
end
%% check boundaries
function [x] = bound(x,lb,ub,len,margin)
% upper bound - build column matrix of bounds
temp = ub * ones(1,len);
% offset from the bound by the margin
junk = (1+margin)*temp - margin*x;
% set values above the upper bound to the limit
x(x > temp) = junk(x > temp);
% lower bound - build column matrix of bounds
temp = lb * ones(1,len);
% offset from the bound by the margin
junk = (1+margin)*temp - margin*x;
% set values below lower bound to the limit
x(x < temp) = junk(x < temp);

function [bus, gen, branch, f, success, et] = opf(baseMVA, bus, gen, gencost, ...
					branch, areas, Ybus, Yf, Yt, ref, pv, pq, mpopt)
%OPF  Solves an optimal power flow.
%   [bus, gen, branch, f, success, et] = opf(baseMVA, bus, gen, gencost, ...
%                   branch, areas, Ybus, Yf, Yt, ref, pv, pq, mpopt)
%   Assumes that Ybus, Yf, Yt, V, ref, pv, pq are consistent with
%   (or override) data in bus, gen, branch.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 1996-2003 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/ for more info.

%%----- initialization -----
%% default arguments
if nargin < 13
	mpopt = mpoption;		%% use default options
end

%% options
verbose	= mpopt(31);
npts = mpopt(14);		%% number of points to evaluate when converting
						%% polynomials to piece-wise linear

%% define constants
j = sqrt(-1);

%% define named indices into data matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
	VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, ...
	GEN_STATUS, PMAX, PMIN, MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN] = idx_gen;
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, N, COST] = idx_cost;

%%-----  check/convert costs, set default algorithm  -----
%% get cost model, check consistency
model = gencost(:, MODEL);
i_pwln = find(model == PW_LINEAR);
i_poly = find(model == POLYNOMIAL);
if any(i_pwln) & any(i_poly) & verbose
	fprintf('not all generators use same cost model, all will be converted to piece-wise linear\n');
end

%% set algorithm
if mpopt(11) == 0
	%% use default for this cost model
	if find(model ~= PW_LINEAR & model ~= POLYNOMIAL)
		error('unknown generator cost model');
	elseif any(i_pwln)		%% some piece-wise linear, use appropriate alg
		mpopt(11) = mpopt(13);
	else					%% must all be polynomial
		mpopt(11) = mpopt(12);
	end
end
alg = mpopt(11);
formulation = opf_form(alg);

% %% move Pmin and Pmax limits out slightly to avoid problems
% %% with lambdas caused by rounding errors when corner point
% %% of cost function lies at exactly Pmin or Pmax
% if any(i_pwln)
% 	ng = size(gen, 1);
% 	gen(:, PMIN) = gen(:, PMIN) - 10 * mpopt(16) * ones(ng, 1);
% 	gen(:, PMAX) = gen(:, PMAX) + 10 * mpopt(16) * ones(ng, 1);
% end

%% check cost model/algorithm consistency
if any( i_pwln ) & formulation == 1
	error(sprintf('algorithm %d does not handle piece-wise linear cost functions', alg));
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%	Eventually put code here to fit polynomials to piece-wise linear as needed.
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%% convert polynomials to piece-wise linear
if any(i_poly)  & formulation == 2
	if verbose
		fprintf('converting from polynomial to piece-wise linear cost model\n');
	end
	[pcost, qcost] = pqcost(gencost, size(gen, 1));
	i_poly = find(pcost(:, MODEL) == POLYNOMIAL);
	pcost = poly2pwl(pcost(i_poly, :), gen(i_poly, PMIN), gen(i_poly, PMAX), npts);
	if ~isempty(qcost)
		i_poly = find(qcost(:, MODEL) == POLYNOMIAL);
		qcost = poly2pwl(qcost(i_poly, :), gen(i_poly, QMIN), gen(i_poly, QMAX), npts);
	end
	gencost = [pcost; qcost];
end

%%-----  run opf  -----
if formulation == 5
	if mpopt(61) == 0
		mpopt(61) = mpopt(16);
	end
	if mpopt(62) == 0
		mpopt(62) = mpopt(16);
	end
	[bus, gen, branch, f, info, g, jac, et] = mopf(baseMVA, ...
				bus, gen, branch, areas, gencost, sparse(0,0), [], [], mpopt);
	success = ~info;
else
	%% start timer
	t0 = clock;
	
	%% gen info
	on = find(gen(:, GEN_STATUS) > 0);		%% which generators are on?
	gbus = gen(on, GEN_BUS);				%% what buses are they at?
	
	%% sizes of things
	nb = size(bus, 1);
	nl = size(branch, 1);
	npv	= length(pv);
	npq	= length(pq);
	ng = length(on);						%% number of generators that are turned on
	
	%% initial state
	V	= bus(:, VM) .* exp(sqrt(-1) * pi/180 * bus(:, VA));
	V(gbus) = gen(on, VG) ./ abs(V(gbus)).* V(gbus);
	Pg	= gen(on, PG) / baseMVA;
	Qg	= gen(on, QG) / baseMVA;
	
	%% check for costs for Qg
	[pcost, qcost] = pqcost(gencost, size(gen, 1), on);
	
	%% set up indexing for x
	j1 = 1;			j2	= npv;				%% j1:j2	- V angle of pv buses
	j3 = j2 + 1;	j4	= j2 + npq;			%% j3:j4	- V angle of pq buses
	j5 = j4 + 1;	j6	= j4 + nb;			%% j5:j6	- V mag of all buses
	j7 = j6 + 1;	j8	= j6 + ng;			%% j7:j8	- P of generators
	j9 = j8 + 1;	j10	= j8 + ng;			%% j9:j10	- Q of generators
	j11 = j10 + 1;	j12	= j10 + ng;			%% j11:j12	- Cp, cost of Pg
	j13 = j12 + 1;	j14	= j12 + ng;			%% j13:j14	- Cq, cost of Qg
	
	%% set up x
	if formulation == 1
		Cp = [];
		Cq = [];
	else
		Cp = totcost(pcost, Pg * baseMVA);
		Cq = totcost(qcost, Qg * baseMVA);	%% empty if qcost is empty
	end
	x = [angle(V([pv; pq])); abs(V); Pg; Qg; Cp; Cq];
	
	%% run constrained optimization
	[fun, grad]	= fg_names(alg);
	mpopt(15)	= 2 * nb;				%% set number of equality constraints
	
	if opf_slvr(alg) == 0			%% use CONSTR
		%% set some options
		if mpopt(19) == 0
			mpopt(19) = 2 * nb + 150;	%% set max number of iterations for constr
		end
	
		%% set up options for Optim Tbx's constr
		otopt = foptions;				%% get default options for constr
		otopt(1) = (verbose > 0);		%% set verbose flag appropriately
		% otopt(9) = 1;					%% check user supplied gradients?
		otopt(2)  = mpopt(17);			%% termination tolerance on 'x'
		otopt(3)  = mpopt(18);			%% termination tolerance on 'F'
		otopt(4)  = mpopt(16);			%% termination tolerance on constraint violation
		otopt(13) = mpopt(15);			%% number of equality constraints
		otopt(14) = mpopt(19);			%% maximum number of iterations
		
		%% run optimization
		[x, otopt, lambda] = constr(fun, x, otopt, [], [], grad, ...
							baseMVA, bus, gen, gencost, branch, Ybus, Yf, Yt, V, ref, pv, pq, mpopt);
		
		%% get final objective function value & constraint values
		[f, g] = feval(fun, x, baseMVA, bus, gen, gencost, branch, Ybus, Yf, Yt, V, ref, pv, pq, mpopt);
		
		%% check for convergence
		if otopt(10) >= otopt(14)	| max(abs(g(1:otopt(13))))			> otopt(4) ...
									| max(g((otopt(13)+1):length(g)))	> otopt(4)
			success = 0;				%% did NOT converge
		else
			success = 1;				%% DID converge
		end
	else							%% use LPCONSTR
		%% run load flow to get starting point
		[x, success_lf] = LPeqslvr(x, baseMVA, bus, gen, gencost, branch, Ybus, Yf, Yt, V, ref, pv, pq, mpopt);
		if success_lf ~= 1
			error('Sorry, cannot find a starting point using power flow, please check data!'); 
		end
		
		%% set step size
		cstep = 0;
		if ~isempty(Cp)
			cstep = max(abs(Cp));
			if cstep < 1.0e6, cstep = 1.0e6; end
		end
		step0=[2*ones(nb-1,1);					%% starting stepsize for Vangle
				ones(nb,1);						%% Vmag
				0.6*ones(ng,1);					%% Pg
				0.3*ones(ng,1);					%% Qg
				cstep*ones(length(Cp),1);		%% Cp
				cstep*ones(length(Cq),1)];		%% Cq
		idx_xi = [];
	
		%% run optimization
		[x, lambda, success] = LPconstr(fun, x, idx_xi, mpopt, step0, [], [], grad, 'LPeqslvr', ...
							baseMVA, bus, gen, gencost, branch, Ybus, Yf, Yt, V, ref, pv, pq, mpopt);
	
		%% get final objective function value & constraint values
		f = feval(fun, x, baseMVA, bus, gen, gencost, branch, Ybus, Yf, Yt, V, ref, pv, pq, mpopt);
	end
	
	%% reconstruct V
	Va = zeros(nb, 1);
	Va([ref; pv; pq]) = [angle(V(ref)); x(j1:j2); x(j3:j4)];
	Vm = x(j5:j6);
	V = Vm .* exp(j * Va);
	
	%% grab Pg & Qg
	Sg = x(j7:j8) + j * x(j9:j10);		%% complex power generation in p.u.
	
	%%-----  calculate return values  -----
	%% update bus, gen, branch with solution info
	[bus, gen, branch] = opfsoln(baseMVA, bus, gen, branch, ...
							Ybus, Yf, Yt, V, Sg, lambda, ref, pv, pq, mpopt);
	
	%% compute elapsed time
	et = etime(clock, t0);
end

return;

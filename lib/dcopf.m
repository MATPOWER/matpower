function [bus, gen, branch, f, success, et] = dcopf(baseMVA, bus, gen, gencost, ...
					branch, B, Bf, Pbusinj, Pfinj, ref, pv, pq, mpopt)
%DCOPF  Solves a DC optimal power flow.
%   [bus, gen, branch, f, success, et] = dcopf(baseMVA, bus, gen, gencost, ...
%                   branch, B, Bf, Pbusinj, Pfinj, ref, pv, pq, mpopt)
%   Assumes that B, Bf, Pbusinj, Pfinj, V, ref, pv, pq are consistent with
%   (or override) data in bus, gen, branch.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   Copyright (c) 1996-2003 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/ for more info.

%%----- initialization -----
%% default arguments
if nargin < 12
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
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, ...
	RATE_C, TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST] = idx_brch;
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
if any(i_pwln) | any(find(gencost(:, N) > 3))
	formulation = 2;	%% use piecewise linear formulation
else
	formulation = 1;	%% use polynomial cost formulation
	if any(find(gencost(:, N) ~= 3))
		error('DC opf with polynomial costs can only handle quadratic costs.');
	end
end

%% move Pmin and Pmax limits out slightly to avoid problems
%% caused by rounding errors when corner point of cost
%% function lies at exactly Pmin or Pmax
if any(i_pwln)
	ng = size(gen, 1);
	gen(:, PMIN) = gen(:, PMIN) - 10 * mpopt(16) * ones(ng, 1);
	gen(:, PMAX) = gen(:, PMAX) + 10 * mpopt(16) * ones(ng, 1);
end

%% convert polynomials to piece-wise linear (if necessary)
if any(i_poly) & formulation == 2
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
%% start timer
t0 = clock;

%% gen info
on = find(gen(:, GEN_STATUS) > 0);		%% which generators are on?
gbus = gen(on, GEN_BUS);				%% what buses are they at?

%% sizes of things
nb = size(bus, 1);
nl = size(branch, 1);
ng = length(on);						%% number of generators that are turned on

%% initial state
Va	= bus(:, VA) * (pi/180);
Pg	= gen(on, PG) / baseMVA;

%% check for costs for Qg
[pcost, qcost] = pqcost(gencost, size(gen, 1), on);

%% set up x along with indexing
j1 = 1;			j2	= nb;				%% j1:j2	- bus V angles
j3 = j2 + 1;	j4	= j2 + ng;			%% j3:j4	- P of generators
j5 = j4 + 1;	j6	= j4 + ng;			%% j5:j6	- Cp, cost of Pg
if formulation == 2				%% piece-wise linear costs
	Cp = totcost(pcost, Pg * baseMVA);
else							%% polynomial costs
	Cp = [];
end
nc = length(Cp);						%% number of cost variables (ng or 0)
x = [Va; Pg; Cp];

%% set up constraints and indexing where,	AA * x <= bb
% i0 = 1;							%% 1 - voltage angle reference
i1 = 2;			i2 = nb + 1;		%% i1:i2 - P mismatch, all buses
i3 = i2 + 1;	i4 = i2 + ng;		%% i3:i4 - Pmin, gen buses
i5 = i4 + 1;	i6 = i4 + ng;		%% i5:i6 - Pmax, gen buses
i7 = i6 + 1;	i8 = i6 + nl;		%% i7:i8 - |Pf| line limit
i9 = i8 + 1;	i10 = i8 + nl;		%% i9:i10 - |Pt| line limit
if formulation == 2				%% piece-wise linear costs
	%% compute cost constraints [ Cp >= m * Pg + b ] => [ m * Pg - Cp <= -b ]
	nsegs = pcost(:, N) - 1;			%% number of cost constraints for each gen
	ncc = sum(nsegs);					%% total number of cost constraints
	Acc = sparse(ncc, nb+ng+nc);
	bcc = zeros(ncc, 1);
	for i = 1:ng
		xx = pcost(i,		COST:2:( COST + 2*(nsegs(i))	))';
		yy = pcost(i,	(COST+1):2:( COST + 2*(nsegs(i)) + 1))';
		k1 = 1:nsegs(i);
		k2 = 2:(nsegs(i) + 1);
		m = (yy(k2) - yy(k1)) ./ (xx(k2) - xx(k1));
		b = yy(k1) - m .* xx(k1);
		Acc(sum(nsegs(1:(i-1))) + [1:nsegs(i)], nb+i)		= m * baseMVA;
		Acc(sum(nsegs(1:(i-1))) + [1:nsegs(i)], nb+ng+i)	= -ones(nsegs(i), 1);
		bcc(sum(nsegs(1:(i-1))) + [1:nsegs(i)])				= -b;
	end
else							%% polynomial costs
	Acc = sparse(0, nb+ng+nc);
	bcc = zeros(0, 1);
end

AA = [
		sparse(1, ref, 1, 1, nb+ng+nc);						%% reference angle
		B,	-sparse(gen(on, GEN_BUS), 1:ng, ones(ng, 1), nb, ng), sparse(nb, nc);
															%% real power flow eqns
		sparse(ng, nb), -speye(ng, ng), sparse(ng, nc);		%% lower limit on Pg
		sparse(ng, nb), speye(ng, ng), sparse(ng, nc);		%% upper limit on Pg
		Bf, sparse(nl, ng+nc);								%% flow limit on Pf
		-Bf, sparse(nl, ng+nc);								%% flow limit on Pt
		Acc;												%% cost constraints
];

bb = [
		Va(ref);											%% reference angle
		-(bus(:, PD) + bus(:, GS)) / baseMVA - Pbusinj;		%% real power flow eqns
		-gen(on, PMIN) / baseMVA;							%% lower limit on Pg
		gen(on, PMAX) / baseMVA;							%% upper limit on Pg
		branch(:, RATE_A) / baseMVA - Pfinj;				%% flow limit on Pf
		branch(:, RATE_A) / baseMVA + Pfinj;				%% flow limit on Pt
		bcc;												%% cost constraints
];

%% set up objective function of the form:  0.5 * x'*H*x + c'*x
if formulation == 2				%% piece-wise linear costs
	H = sparse(nb+ng+nc, nb+ng+nc);
	c = [	zeros(nb+ng, 1);  
			ones(nc, 1)		];
else							%% polynomial costs
	polycf = pcost(:, COST:COST+2) * diag([ baseMVA^2 baseMVA 1]);  %% coeffs for Pg in p.u.
	H = sparse(j3:j4, j3:j4, 2*polycf(:, 1), nb+ng, nb+ng );
	c = [	zeros(nb, 1);  
			polycf(:, 2)	];
end

%% run QP solver
mpopt(15)	= nb + 1;				%% set number of equality constraints
if verbose > 1						%% print QP progress for verbose levels 2 & 3
	qpverbose = 1;
else
	qpverbose = -1;
end
if mpopt(51) == 0					%% QP solver can't handle sparse matrices?
	AA = full(AA);
	H = full(H);
end

[x, lambda, how] = qp(H, c, AA, bb, [], [], x, mpopt(15), qpverbose);

if how(1) ~= 'o'
	success = 0;
else
	success = 1;
end

%% update solution data
Va = x(j1:j2);
Pg = x(j3:j4);
f = sum(totcost(pcost, Pg * baseMVA));

%%-----  calculate return values  -----
%% update voltages & generator outputs
bus(:, VM) = ones(nb, 1);
bus(:, VA) = Va * 180 / pi;
gen(on, PG) = Pg * baseMVA;

%% compute flows etc.
branch(:, [QF, QT]) = zeros(nl, 2);
branch(:, PF) = (Bf * Va + Pfinj) * baseMVA;
branch(:, PT) = -branch(:, PF);

%% update lambda's and mu's
bus(:, [LAM_P, LAM_Q, MU_VMIN, MU_VMAX]) = zeros(nb, 4);
gen(:, [MU_PMIN, MU_PMAX, MU_QMIN, MU_QMAX]) = zeros(size(gen, 1), 4);
branch(:, [MU_SF, MU_ST]) = zeros(nl, 2);
bus(:, LAM_P)		= lambda(i1:i2) / baseMVA;
gen(on, MU_PMIN)	= lambda(i3:i4) / baseMVA;
gen(on, MU_PMAX)	= lambda(i5:i6) / baseMVA;
branch(:, MU_SF)	= lambda(i7:i8) / baseMVA;
branch(:, MU_ST)	= lambda(i9:i10) / baseMVA;

%% compute elapsed time
et = etime(clock, t0);

return;

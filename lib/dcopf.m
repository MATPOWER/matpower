function [buso, geno, brancho, f, success, info, et] = dcopf(baseMVA, ...
                   bus, gen, branch, areas, gencost, mpopt)
%DCOPF  Solves a DC optimal power flow.
%
%   [bus, gen, branch, f, success] = dcopf(casefile, mpopt)
%
%   [bus, gen, branch, f, success] = dcopf(baseMVA, bus, gen, branch, areas, ...
%                                    gencost, mpopt)
%
%   [bus, gen, branch, f, success, info, et] = dcopf(casefile)
%
%   The data for the problem can be specified in one of 3 ways: (1) the name of
%   a case file which defines the data matrices baseMVA, bus, gen, branch,
%   areas and gencost, (2) a struct containing the data matrices as fields, or
%   (3) the data matrices themselves.
%
%   The optional mpopt vector specifies MATPOWER options. Type 'help mpoption'
%   for details and default values.
%
%   The solved case is returned in the data matrices, bus, gen and branch. Also,
%   returned are the final objective function value (f) and a flag which is
%   true if the algorithm was successful in finding a solution (success).
%   Additional optional return values are an algorithm specific return status
%   (info), elapsed time in seconds (et).

%   MATPOWER
%   $Id$
%   by Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Autonoma de Manizales
%   and Ray Zimmerman, PSERC Cornell
%   Copyright (c) 1996-2005 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%%----- initialization -----
%% Sort arguments
if isstr(baseMVA) | isstruct(baseMVA)
  casefile = baseMVA;
  if nargin == 1
    mpopt = mpoption;
  elseif nargin == 2
    mpopt = bus;
  else
    error('dcopf.m: Incorrect input parameter order, number or type');
  end
  [baseMVA, bus, gen, branch, areas, gencost] = loadcase(casefile);
else
  if nargin == 6
    mpopt = mpoption;
  elseif nargin ~= 7
    error('dcopf.m: Incorrect input parameter order, number or type');
  end
end

%% options
mpopt(10) = 1;  % force DC treatment
verbose = mpopt(31);
npts = mpopt(14);       %% number of points to evaluate when converting
                        %% polynomials to piece-wise linear

%% define constants
j = sqrt(-1);

%% define named indices into data matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, ...
    PMAX, PMIN, MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, QMAX2, QMIN2, ...
    RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q] = idx_gen;
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, N, COST] = idx_cost;

% If tables do not have multiplier/extra columns, append zero cols
if size(bus,2) < MU_VMIN
  bus = [bus zeros(size(bus,1),MU_VMIN-size(bus,2)) ];
end
if size(gen,2) < MU_QMIN
  gen = [ gen zeros(size(gen,1),MU_QMIN-size(gen,2)) ];
end
if size(branch,2) < MU_ST
  branch = [ branch zeros(size(branch,1),MU_ST-size(branch,2)) ];
end

% Filter out inactive generators and branches; save original bus & branch
comgen = find(gen(:,GEN_STATUS) > 0);
offgen = find(gen(:,GEN_STATUS) <= 0);
onbranch  = find(branch(:,BR_STATUS) ~= 0);
offbranch = find(branch(:,BR_STATUS) == 0);
genorg = gen;
branchorg = branch;
ng = size(gen,1);         % original size(gen), at least temporally
gen   = gen(comgen, :);
branch = branch(onbranch, :);
if size(gencost,1) == ng
  gencost = gencost(comgen, :);
else
  gencost = gencost( [comgen; comgen+ng], :);
end

% Renumber buses consecutively
[i2e, bus, gen, branch, areas] = ext2int(bus, gen, branch, areas);

%% get bus index lists of each type of bus
[ref, pv, pq] = bustypes(bus, gen);

%% build B matrices and phase shift injections
[B, Bf, Pbusinj, Pfinj] = makeBdc(baseMVA, bus, branch);

%%-----  check/convert costs, set default algorithm  -----
%% get cost model, check consistency
model = gencost(:, MODEL);
i_pwln = find(model == PW_LINEAR);
i_poly = find(model == POLYNOMIAL);
if any(i_pwln) & any(i_poly) & verbose
    fprintf('not all generators use same cost model, all will be converted to piece-wise linear\n');
end
if any(i_pwln) | any(find(gencost(:, N) > 3))
    formulation = 2;    %% use piecewise linear formulation
else
    formulation = 1;    %% use polynomial cost formulation
    if any(find(gencost(:, N) ~= 3))
        error('DC opf with polynomial costs can only handle quadratic costs.');
    end
end

%% convert polynomials to piece-wise linear (if necessary)
if any(i_poly) & formulation == 2
    if verbose
        fprintf('converting from polynomial to piece-wise linear cost model\n');
    end
    [pcost, qcost] = pqcost(gencost, size(gen, 1));
    i_poly = find(pcost(:, MODEL) == POLYNOMIAL);
    tmp = poly2pwl(pcost(i_poly, :), gen(i_poly, PMIN), gen(i_poly, PMAX), npts);
    pcost(i_poly, 1:size(tmp,2)) = tmp;
    if ~isempty(qcost)
        i_poly = find(qcost(:, MODEL) == POLYNOMIAL);
        tmp = poly2pwl(qcost(i_poly, :), gen(i_poly, QMIN), gen(i_poly, QMAX), npts);
        qcost(i_poly, 1:size(tmp,2)) = tmp;
    end
    gencost = [pcost; qcost];
end

%%-----  run opf  -----
%% start timer
t0 = clock;

%% gen info
on = find(gen(:, GEN_STATUS) > 0);      %% which generators are on?
gbus = gen(on, GEN_BUS);                %% what buses are they at?

%% sizes of things
nb = size(bus, 1);
nl = size(branch, 1);
ng = length(on);                        %% number of generators that are turned on

%% initial state
Va  = bus(:, VA) * (pi/180);
Pg  = gen(on, PG) / baseMVA;

%% check for costs for Qg
[pcost, qcost] = pqcost(gencost, size(gen, 1), on);

%% set up x along with indexing
j1 = 1;         j2  = nb;               %% j1:j2    - bus V angles
j3 = j2 + 1;    j4  = j2 + ng;          %% j3:j4    - P of generators
j5 = j4 + 1;    j6  = j4 + ng;          %% j5:j6    - Cp, cost of Pg
if formulation == 2             %% piece-wise linear costs
    Cp = totcost(pcost, Pg * baseMVA);
else                            %% polynomial costs
    Cp = [];
end
nc = length(Cp);                        %% number of cost variables (ng or 0)
x = [Va; Pg; Cp];

%% set up constraints and indexing where,   AA * x <= bb
% i0 = 1;                           %% 1 - voltage angle reference
i1 = 2;         i2 = nb + 1;        %% i1:i2 - P mismatch, all buses
i3 = i2 + 1;    i4 = i2 + ng;       %% i3:i4 - Pmin, gen buses
i5 = i4 + 1;    i6 = i4 + ng;       %% i5:i6 - Pmax, gen buses
i7 = i6 + 1;    i8 = i6 + nl;       %% i7:i8 - |Pf| line limit
i9 = i8 + 1;    i10 = i8 + nl;      %% i9:i10 - |Pt| line limit
if formulation == 2             %% piece-wise linear costs
    %% compute cost constraints [ Cp >= m * Pg + b ] => [ m * Pg - Cp <= -b ]
    nsegs = pcost(:, N) - 1;            %% number of cost constraints for each gen
    ncc = sum(nsegs);                   %% total number of cost constraints
    Acc = sparse(ncc, nb+ng+nc);
    bcc = zeros(ncc, 1);
    for i = 1:ng
        xx = pcost(i,       COST:2:( COST + 2*(nsegs(i))    ))';
        yy = pcost(i,   (COST+1):2:( COST + 2*(nsegs(i)) + 1))';
        k1 = 1:nsegs(i);
        k2 = 2:(nsegs(i) + 1);
        m = (yy(k2) - yy(k1)) ./ (xx(k2) - xx(k1));
        b = yy(k1) - m .* xx(k1);
        Acc(sum(nsegs(1:(i-1))) + [1:nsegs(i)], nb+i)       = m * baseMVA;
        Acc(sum(nsegs(1:(i-1))) + [1:nsegs(i)], nb+ng+i)    = -ones(nsegs(i), 1);
        bcc(sum(nsegs(1:(i-1))) + [1:nsegs(i)])             = -b;
    end
else                            %% polynomial costs
    Acc = sparse(0, nb+ng+nc);
    bcc = zeros(0, 1);
end

AA = [
        sparse(1, ref, 1, 1, nb+ng+nc);                     %% reference angle
        B,  -sparse(gen(on, GEN_BUS), 1:ng, ones(ng, 1), nb, ng), sparse(nb, nc);
                                                            %% real power flow eqns
        sparse(ng, nb), -speye(ng, ng), sparse(ng, nc);     %% lower limit on Pg
        sparse(ng, nb), speye(ng, ng), sparse(ng, nc);      %% upper limit on Pg
        Bf, sparse(nl, ng+nc);                              %% flow limit on Pf
        -Bf, sparse(nl, ng+nc);                             %% flow limit on Pt
        Acc;                                                %% cost constraints
];

bb = [
        Va(ref);                                            %% reference angle
        -(bus(:, PD) + bus(:, GS)) / baseMVA - Pbusinj;     %% real power flow eqns
        -gen(on, PMIN) / baseMVA;                           %% lower limit on Pg
        gen(on, PMAX) / baseMVA;                            %% upper limit on Pg
        branch(:, RATE_A) / baseMVA - Pfinj;                %% flow limit on Pf
        branch(:, RATE_A) / baseMVA + Pfinj;                %% flow limit on Pt
        bcc;                                                %% cost constraints
];

%% set up objective function of the form:  0.5 * x'*H*x + c'*x
if formulation == 2             %% piece-wise linear costs
    H = sparse(nb+ng+nc, nb+ng+nc);
    c = [   zeros(nb+ng, 1);  
            ones(nc, 1)     ];
else                            %% polynomial costs
    polycf = pcost(:, COST:COST+2) * diag([ baseMVA^2 baseMVA 1]);  %% coeffs for Pg in p.u.
    H = sparse(j3:j4, j3:j4, 2*polycf(:, 1), nb+ng, nb+ng );
    c = [   zeros(nb, 1);  
            polycf(:, 2)    ];
end

%% run QP solver
mpopt(15)   = nb + 1;               %% set number of equality constraints
if verbose > 1                      %% print QP progress for verbose levels 2 & 3
    qpverbose = 1;
else
    qpverbose = -1;
end
if ~have_fcn('sparse_qp') | mpopt(51) == 0 %% don't use sparse matrices
    AA = full(AA);
    H = full(H);
end

[x, lambda, how, success] = mp_qp(H, c, AA, bb, [], [], x, mpopt(15), qpverbose);

info = success;

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
bus(:, LAM_P)       = lambda(i1:i2) / baseMVA;
gen(on, MU_PMIN)    = lambda(i3:i4) / baseMVA;
gen(on, MU_PMAX)    = lambda(i5:i6) / baseMVA;
branch(:, MU_SF)    = lambda(i7:i8) / baseMVA;
branch(:, MU_ST)    = lambda(i9:i10) / baseMVA;

% convert to original external bus ordering
[bus, gen, branch, areas] = int2ext(i2e, bus, gen, branch, areas);

% Now create output matrices with all lines, all generators, committed and
% non-committed
geno = genorg;
brancho = branchorg;
geno(comgen, : ) = gen;
brancho(onbranch, :)  = branch;
% And zero out appropriate fields of non-comitted generators and lines
tmp = zeros(length(offgen), 1);
geno(offgen, PG) = tmp;
geno(offgen, QG) = tmp;
geno(offgen, MU_PMAX) = tmp;
geno(offgen, MU_PMIN) = tmp;
tmp = zeros(length(offbranch), 1);
brancho(offbranch, PF) = tmp;
brancho(offbranch, QF) = tmp;
brancho(offbranch, PT) = tmp;
brancho(offbranch, QT) = tmp;
brancho(offbranch, MU_SF) = tmp;
brancho(offbranch, MU_ST) = tmp;

%% compute elapsed time
et = etime(clock, t0);

if (nargout == 0) & (info > 0)
  printpf(baseMVA, bus, geno, brancho, f, success, et, 1, mpopt);
end

if nargout, buso = bus; end

return;

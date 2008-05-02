function [busout, genout, branchout, f, success, info, et, xz, pimul] = ...
    dcopf(varargin)
%DCOPF  Solves a DC optimal power flow.
%
%   [bus, gen, branch, f, success] = dcopf(casefile, mpopt)
%
%   [bus, gen, branch, f, success] = dcopf(casefile, A, l, u, mpopt)
%
%   [bus, gen, branch, f, success] = dcopf(baseMVA, bus, gen, branch, ...
%                                    areas, gencost, mpopt)
%
%   [bus, gen, branch, f, success] = dcopf(baseMVA, bus, gen, branch, ...
%                                    areas, gencost, A, l, u, mpopt)
%
%   [bus, gen, branch, f, success] = dcopf(baseMVA, bus, gen, branch, ...
%                                    areas, gencost, A, l, u, mpopt, ...
%                                    N, fparm, H, Cw)
%
%   [bus, gen, branch, f, success] = dcopf(baseMVA, bus, gen, branch, ...
%                                    areas, gencost, A, l, u, mpopt, ...
%                                    N, fparm, H, Cw, z0, zl, zu)
%
%   [bus, gen, branch, f, success, info, et, xr, pimul] = dcopf(casefile)
%
%
%   The data for the problem can be specified in one of 3 ways: (1) the name of
%   a case file which defines the data matrices baseMVA, bus, gen, branch,
%   areas and gencost, (2) a struct containing the data matrices as fields, or
%   (3) the data matrices themselves.
%
%   When specified, A, l, u represent additional linear constraints on the
%   optimization variables, l <= A*[x; z] <= u. For an explanation of the
%   formulation used and instructions for forming the A matrix, type
%   'help genform'.
%
%   A generalized cost on all variables can be applied if input arguments
%   N, fparm, H and Cw are specified.  First, a linear transformation
%   of the optimization variables is defined by means of r = N * [x; z].
%   Then, to each element of r a function is applied as encoded in the
%   fparm matrix (see manual).  If the resulting vector is now named w,
%   then H and Cw define a quadratic cost on w: (1/2)*w'*H*w + Cw * w .
%   H and N should be sparse matrices and H should also be symmetric.
%
%   The optional mpopt vector specifies MATPOWER options. Type 'help mpoption'
%   for details and default values.
%
%   The solved case is returned in the data matrices, bus, gen and branch. Also
%   returned are the final objective function value (f) and a flag which is
%   true if the algorithm was successful in finding a solution (success).
%   Additional optional return values are an algorithm specific return status
%   (info), elapsed time in seconds (et), the vector of variables (xr) as well 
%   as the constraint multipliers (pimul).
%
%   Rules for A matrix: If the user specifies an A matrix that has more columns
%   than the number of "x" (OPF) variables, then there are extra linearly
%   constrained "z" variables.

%   MATPOWER
%   $Id$
%   by Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Autonoma de Manizales
%   and Ray Zimmerman, PSERC Cornell
%   Copyright (c) 1996-2008 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

t0 = clock;

% process input arguments
[baseMVA, bus, gen, branch, areas, gencost, Au, lbu, ubu, mpopt, ...
    N, fparm, H, Cw, z0, zl, zu] = opf_args(varargin{:});

%%----- initialization -----
verbose = mpopt(31);
mpopt(10) = 1;          %% force DC treatment

%% define named indices into data matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;

% If tables do not have multiplier/extra columns, append zero cols.
% Update whenever the data format changes!
if size(bus,2) < MU_VMIN
  bus = [bus zeros(size(bus,1),MU_VMIN-size(bus,2)) ];
end
if size(gen,2) < MU_QMIN
  gen = [ gen zeros(size(gen,1),MU_QMIN-size(gen,2)) ];
end
if size(branch,2) < MU_ANGMAX
  branch = [ branch zeros(size(branch,1),MU_ANGMAX-size(branch,2)) ];
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

%% ignore reactive costs
[gencost, qcost] = pqcost(gencost, size(gen, 1));

% Renumber buses consecutively
[i2e, bus, gen, branch, areas] = ext2int(bus, gen, branch, areas);

%% get bus index lists of each type of bus
[ref, pv, pq] = bustypes(bus, gen);

% Print a warning if there is more than one reference bus
refs = find(bus(:, BUS_TYPE) == REF);
if length(refs) > 1
  errstr = ['\ndcopf: Warning: more than one reference bus detected in bus table data.\n', ...
              '       For a system with islands, a reference bus in each island\n', ...
              '       might help convergence but in a fully connected system such\n', ...
              '       a situation is probably not reasonable.\n\n' ];
  fprintf(errstr);
end

% Find out which branches require angle constraints
if mpopt(25)        % OPF_IGNORE_ANG_LIM
  nang = 0;
else
  iang = find((branch(:, ANGMIN) & branch(:, ANGMIN) > -360) | ...
              (branch(:, ANGMAX) & branch(:, ANGMAX) < 360));
  iangl = find(branch(iang, ANGMIN));
  iangh = find(branch(iang, ANGMAX));
  nang = length(iang);
end

%% problem dimensions
ipol = find(gencost(:, MODEL) == POLYNOMIAL); %% polynomial costs
ipwl = find(gencost(:, MODEL) == PW_LINEAR);  %% piece-wise linear costs
nb = size(bus, 1);      %% number of buses
nl = size(branch, 1);   %% number of branches
ng = size(gen, 1);      %% number of dispatchable injections
nx = nb + ng;           %% number of standard OPF control variables
ny = size(ipwl, 1);     %% number of piece-wise linear costs
nusr = size(Au, 1);     %% number of linear user constraints
nw = size(N, 1);        %% number of general cost vars, w
if isempty(Au)
  nz = 0;
  Au = sparse(0,nx);
  if ~isempty(N)        %% still need to check number of columns of N
    if size(N, 2) ~= nx;
      error(sprintf('dcopf.m: user supplied N matrix must have %d columns.', nx));
    end
  end
else
  nz = size(Au,2) - nx; %% additional linear variables
  if nz < 0
    error(sprintf('dcopf.m: user supplied A matrix must have at least %d columns.', nx));
  end
end
nxyz = nx+ny+nz;        %% total number of control vars of all types

%% define indexing of optimization variable vector
k = 0;
thbas   = k + 1;    k = k + nb;     thend = k;      %% voltage angles
pgbas   = k + 1;    k = k + ng;     pgend = k;      %% real power injections
ybas    = k + 1;    k = k + ny;     yend  = k;      %% pwl costs
zbas    = k + 1;    k = k + nz;     zend  = k;      %% user vars

%% define indexing of constraint vector
k = 0;
pmsmbas = k + 1;    k = k + nb;     pmsmend = k;    %% real power balance
sfbas   = k + 1;    k = k + nl;     sfend   = k;    %% "from" flow limit
stbas   = k + 1;    k = k + nl;     stend   = k;    %% "to" flow limit
usrbas  = k + 1;    k = k + nusr;   usrend  = k;    % warning: nusr could be 0
angbas  = k + 1;    k = k + nang;   angend  = k;    %% branch angle diff lims
%% not done yet, need number of Ay constraints.

% Let makeAy deal with any y-variable for piecewise-linear convex costs.
% note that if there are z variables then Ay doesn't have the columns
% that would span the z variables, so we append them.
if ny > 0
  [Ay, by]  = makeAy(baseMVA, ng, gencost, pgbas, [], ybas);
  if nz > 0
    Ay = [ Ay  sparse(size(Ay,1), nz) ];
  end
else
  Ay = [];
  by =[];
end
ncony = size(Ay,1);
yconbas = k + 1;    k = k + ncony;  yconend = k;    %% done w/constraint idxing

% Make Aang, lang, uang for branch angle difference limits
if nang > 0
  ii = [(1:nang)'; (1:nang)'];
  jj = [branch(iang, F_BUS); branch(iang, T_BUS)];
  Aang = sparse(ii, jj, [ones(nang, 1); -ones(nang, 1)], nang, nxyz);
  uang = 1e10*ones(nang,1);
  lang = -uang;
  lang(iangl) = branch(iang(iangl), ANGMIN) * pi/180;
  uang(iangh) = branch(iang(iangh), ANGMAX) * pi/180;
else
  Aang =[];
  lang =[];
  uang =[];
end

% Insert y columns in Au and N as necessary
if ny > 0
  if nz > 0
    Au = [ Au(:,1:pgend) sparse(nusr, ny) Au(:, pgend+(1:nz)) ];
    if ~isempty(N)
        N = [ N(:,1:pgend) sparse(nw, ny) N(:, pgend+(1:nz)) ];
    end
  else
    Au = [ Au sparse(nusr, ny) ];
    if ~isempty(N)
        N = [ N sparse(nw, ny) ];
    end
  end
end

%% build B matrices and phase shift injections
[B, Bf, Pbusinj, Pfinj] = makeBdc(baseMVA, bus, branch);
neg_Cg = sparse(gen(:, GEN_BUS), 1:ng, -1, nb, ng);   %% Pbus w.r.t. Pg

%% set up constraints as l < A*x < u
A = [ ...
    B,  neg_Cg, sparse(nb, ny+nz);      %% real power balance
    Bf,  sparse(nl, ng+ny+nz);          %% flow limit on Pf
    -Bf, sparse(nl, ng+ny+nz);          %% flow limit on Pt (= -Pf)
    Au;
    Aang;
    Ay
];

l = [ ...
    -(bus(:, PD) + bus(:, GS)) / baseMVA - Pbusinj;
    -1e10 * ones(nl, 1);
    -1e10 * ones(nl, 1);
    lbu;
    lang;
    -1e10 * ones(ncony, 1);
];

u = [ ...
    -(bus(:, PD) + bus(:, GS)) / baseMVA - Pbusinj;
    branch(:, RATE_A) / baseMVA - Pfinj;
    branch(:, RATE_A) / baseMVA + Pfinj;
    ubu;
    uang;
    by;
];

%% set up bounds on optimization vars; y vars unbounded
Va   = bus(:, VA) * (pi/180);
Pg   = gen(:, PG) / baseMVA;
Pmin = gen(:, PMIN) / baseMVA;
Pmax = gen(:, PMAX) / baseMVA;
UB = Inf * ones(nxyz, 1);
LB = -UB;
LB(thbas-1+refs) = Va(refs);        UB(thbas-1+refs) = Va(refs);
LB(pgbas:pgend)  = Pmin;            UB(pgbas:pgend) = Pmax;
if ~isempty(zl)
  LB(zbas:zend) = zl;
end
if ~isempty(zu)
  UB(zbas:zend) = zu;
end

%% initialize optimization vars
x0 = zeros(nxyz, 1);
x0(thbas:thend) = Va;
x0(pgbas:pgend) = Pg;
if ny
  x0(ybas:yend) = totcost(gencost(ipwl, :), Pg * baseMVA);
end
if isempty(z0)
  z0 = zeros(nz, 1);
end
if ~isempty(z0)
  x0(zbas:zend) = z0;
end

%% divide l <= A*x <= u into less than, equal to, greater than,
%% doubly-bounded sets
ieq = find( abs(u-l) <= eps);       %% equality
igt = find( u >= 1e10);             %% greater than, unbounded above
ilt = find( l <= -1e10);            %% less than, unbounded below
ibx = find( (abs(u-l) > eps) & (u < 1e10) & (l > -1e10));
% Aeq = A(ieq, :);
% beq = u(ieq);
% Ai  = [ A(ilt, :); -A(igt, :); A(ibx, :); -A(ibx, :) ];
% bi  = [ u(ilt);   -l(igt);     u(ibx);    -l(ibx)];
AA  = [ A(ieq, :);  A(ilt, :);  -A(igt, :);  A(ibx, :);  -A(ibx, :) ];
bb  = [ u(ieq);     u(ilt);     -l(igt);     u(ibx);     -l(ibx)    ];

%% set up objective function of the form: f = 1/2 * X'*HH*X + CC'*X
%% where X = [x;y;z]. First set up as quadratic function of w,
%% f = 1/2 * w'*HHw*w + CCw'*w, where w = diag(M) * (N*X - Rhat). We
%% will be building on the (optionally present) user supplied parameters.

%% piece-wise linear costs
npwl = length(ipwl);
Npwl = sparse(ones(npwl,1), ybas-1+ipwl, 1, 1, nxyz);       %% sum of y vars
Hpwl = 0;
Cpwl = 1;
fparm_pwl = [1 0 0 1];

%% quadratic costs
npol = length(ipol);
if any(find(gencost(ipol, NCOST) > 3))
    error('DC opf cannot handle polynomial costs with higher than quadratic order.');
end
iqdr = find(gencost(ipol, NCOST) == 3);
ilin = find(gencost(ipol, NCOST) == 2);
polycf = zeros(npol, 3);                            %% quadratic coeffs for Pg
polycf(iqdr, :)   = gencost(ipol(iqdr), COST:COST+2);
polycf(ilin, 2:3) = gencost(ipol(ilin), COST:COST+1);
polycf = polycf * diag([ baseMVA^2 baseMVA 1]);     %% convert to p.u.
Npol = sparse(1:npol, pgbas-1+ipol, 1, npol, nxyz);         %% Pg vars
Hpol = sparse(1:npol, 1:npol, 2*polycf(:, 1), npol, npol);
Cpol = polycf(:, 2);
fparm_pol = ones(npol,1) * [ 1 0 0 1 ];

%% combine with user costs
NN = [ Npwl; Npol; N ];
HHw = [ Hpwl, sparse(1, npol+nw);
        sparse(npol, 1), Hpol, sparse(npol, nw);
        sparse(nw, 1+npol), H   ];
CCw = [Cpwl; Cpol; Cw];
ffparm = [ fparm_pwl; fparm_pol; fparm ];

%% transform quadratic coefficients for w into coefficients for X
nnw = 1+npol+nw;
M   = spdiags(ffparm(:, 4), 0, nnw, nnw);
MR  = M * ffparm(:, 2);
HMR = HHw * MR;
MN  = M * NN;
HH = MN' * HHw * MN;
CC = full(MN' * (CCw - HMR));
C0 = 1/2 * MR' * HMR + sum(polycf(:, 3));   %% constant term of cost

%% run QP solver
mpopt(15) = length(ieq);            %% set number of equality constraints
if verbose > 1                      %% print QP progress for verbose levels 2 & 3
    qpverbose = 1;
else
    qpverbose = -1;
end
if ~have_fcn('sparse_qp') | mpopt(51) == 0 %% don't use sparse matrices
    AA = full(AA);
    HH = full(HH);
end

%%-----  run opf  -----
if any(any(HH))
  [x, lambda, how, success] = mp_qp(HH, CC, AA, bb, LB, UB, x0, mpopt(15), qpverbose, 0);
else
  [x, lambda, how, success] = mp_lp(CC, AA, bb, LB, UB, x0, mpopt(15), qpverbose, 0);
end
info = success;

%% update solution data
Va = x(thbas:thend);
Pg = x(pgbas:pgend);
z  = x(zbas:zend);
f = 1/2 * x' * HH * x + CC' * x + C0;

%%-----  calculate return values  -----
%% update voltages & generator outputs
bus(:, VM) = ones(nb, 1);
bus(:, VA) = Va * 180 / pi;
gen(:, PG) = Pg * baseMVA;

%% compute branch flows
branch(:, [QF, QT]) = zeros(nl, 2);
branch(:, PF) = (Bf * Va + Pfinj) * baseMVA;
branch(:, PT) = -branch(:, PF);

% extract lambdas from constraints
nA = length(u);
neq = length(ieq);
nlt = length(ilt);
ngt = length(igt);
nbx = length(ibx);
lam = zeros(nA, 1);
lam(ieq) = lambda(1:neq);
lam(ilt) = lambda(neq+[1:nlt]);
lam(igt) = -lambda(neq+nlt+[1:ngt]);
lam(ibx) = lambda(neq+nlt+ngt+[1:nbx]) - lambda(neq+nlt+ngt+nbx+[1:nbx]);
muLB = lambda(nA+[1:nxyz]);
muUB = lambda(nA+nxyz+[1:nxyz]);

%% update lambda's and mu's
bus(:, [LAM_P, LAM_Q, MU_VMIN, MU_VMAX]) = zeros(nb, 4);
gen(:, [MU_PMIN, MU_PMAX, MU_QMIN, MU_QMAX]) = zeros(size(gen, 1), 4);
branch(:, [MU_SF, MU_ST]) = zeros(nl, 2);
bus(:, LAM_P)       = lam(pmsmbas:pmsmend) / baseMVA;
branch(:, MU_SF)    = lam(sfbas:sfend) / baseMVA;
branch(:, MU_ST)    = lam(stbas:stend) / baseMVA;
gen(:, MU_PMIN)     = muLB(pgbas:pgend) / baseMVA;
gen(:, MU_PMAX)     = muUB(pgbas:pgend) / baseMVA;

% convert to original external bus ordering
[bus, gen, branch, areas] = int2ext(i2e, bus, gen, branch, areas);

% stick in non-linear constraints too, so we can use the indexing variables
% we've defined, and negate so it looks like the pimul from MINOS
pimul = [
  -lam(pmsmbas:pmsmend);
  -lam(sfbas:sfend);
  -lam(stbas:stend);
  -lam(usrbas:yconend);
  -1;       %% dummy entry corresponding to linear cost row in A (in MINOS)
  muLB - muUB
];

% Now create output matrices with all lines, all generators, committed and
% non-committed
genout = genorg;
branchout = branchorg;
genout(comgen, : ) = gen;
branchout(onbranch, :)  = branch;
% And zero out appropriate fields of non-comitted generators and lines
if ~isempty(offgen)
  tmp = zeros(length(offgen), 1);
  genout(offgen, PG) = tmp;
  genout(offgen, QG) = tmp;
  genout(offgen, MU_PMAX) = tmp;
  genout(offgen, MU_PMIN) = tmp;
end
if ~isempty(offbranch)
  tmp = zeros(length(offbranch), 1);
  branchout(offbranch, PF) = tmp;
  branchout(offbranch, QF) = tmp;
  branchout(offbranch, PT) = tmp;
  branchout(offbranch, QT) = tmp;
  branchout(offbranch, MU_SF) = tmp;
  branchout(offbranch, MU_ST) = tmp;
end

%% compute elapsed time
et = etime(clock, t0);

if (nargout == 0) & (info > 0)
  printpf(baseMVA, bus, genout, branchout, f, success, et, 1, mpopt);
end

if nargout, busout = bus; end

return;

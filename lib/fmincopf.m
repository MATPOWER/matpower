function [busout, genout, branchout, f, success, info, et, g, jac, x, pimul] = ...
      fmincopf(varargin)
%FMINCOPF  Solves an AC optimal power flow using FMINCON (Opt Tbx 2.x & later).
%
%   [bus, gen, branch, f, success] = fmincopf(casefile, mpopt)
%
%   [bus, gen, branch, f, success] = fmincopf(casefile, A, l, u, mpopt)
%
%   [bus, gen, branch, f, success] = fmincopf(baseMVA, bus, gen, branch, ...
%                                    areas, gencost, mpopt)
%
%   [bus, gen, branch, f, success] = fmincopf(baseMVA, bus, gen, branch, ...
%                                    areas, gencost, A, l, u, mpopt)
%
%   [bus, gen, branch, f, success] = fmincopf(baseMVA, bus, gen, branch, ...
%                                    areas, gencost, A, l, u, mpopt, ...
%                                    N, fparm, H, Cw)
%
%   [bus, gen, branch, f, success] = fmincopf(baseMVA, bus, gen, branch, ...
%                                    areas, gencost, A, l, u, mpopt, ...
%                                    N, fparm, H, Cw, z0, zl, zu)
%
%   [bus, gen, branch, f, success, info, et, g, jac, xr, pimul] = fmincopf(casefile)
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
%   (info), elapsed time in seconds (et), the constraint vector (g), the
%   Jacobian matrix (jac), and the vector of variables (xr) as well 
%   as the constraint multipliers (pimul).
%
%   Rules for A matrix: If the user specifies an A matrix that has more columns
%   than the number of "x" (OPF) variables, then there are extra linearly
%   constrained "z" variables.
%
%   NOTE: The shadow prices (lambda's and mu's) produced by some versions of
%         fmincon appear to be slightly inaccurate.

%   MATPOWER
%   $Id$
%   by Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Autonoma de Manizales
%   Copyright (c) 2000-2006 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%%----- initialization -----
t0 = clock;         %% start timer
% process input arguments
[baseMVA, bus, gen, branch, areas, gencost, Au, lbu, ubu, mpopt, ...
    N, fparm, H, Cw, z0, zl, zu] = opf_args(varargin{:});

%% options
verbose = mpopt(31);    %% VERBOSE

%% define constants
j = sqrt(-1);

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

% Renumber buses consecutively
[i2e, bus, gen, branch, areas] = ext2int(bus, gen, branch, areas);

% Sort generators in order of increasing bus number;
ng = size(gen,1);
[tmp, igen] = sort(gen(:, GEN_BUS));
[tmp, inv_gen_ord] = sort(igen);  % save for inverse reordering at the end
gen  = gen(igen, :);
if ng == size(gencost,1)
  gencost = gencost(igen, :);
else
  gencost = gencost( [igen; igen+ng], :);
end

% Print a warning if there is more than one reference bus
refs = find(bus(:, BUS_TYPE) == REF);
if length(refs) > 1
  errstr = ['\nfmincopf: Warning: more than one reference bus detected in bus table data.\n', ...
              '          For a system with islands, a reference bus in each island\n', ...
              '          might help convergence but in a fully connected system such\n', ...
              '          a situation is probably not reasonable.\n\n' ];
  fprintf(errstr);
end

% Find out if any of these "generators" are actually dispatchable loads.
% (see 'help isload' for details on what constitutes a dispatchable load)
% Dispatchable loads are modeled as generators with an added constant
% power factor constraint. The power factor is derived from the original
% value of Pmin and either Qmin (for inductive loads) or Qmax (for capacitive
% loads). If both Qmin and Qmax are zero, this implies a unity power factor
% without the need for an additional constraint.
vload = find( isload(gen) & (gen(:, QMIN) ~= 0 | gen(:, QMAX) ~= 0) );
% At least one of the Q limits must be zero (corresponding to Pmax == 0)
if any( gen(vload, QMIN) ~= 0 & gen(vload, QMAX) ~= 0 )
    error('fmincopf.m: Either Qmin or Qmax must be equal to zero for each dispatchable load.');
end
% Initial values of PG and QG must be consistent with specified power factor
% This is to prevent a user from unknowingly using a case file which would
% have defined a different power factor constraint under a previous version
% which used PG and QG to define the power factor.
Qlim = (gen(vload, QMIN) == 0) .* gen(vload, QMAX) + ...
    (gen(vload, QMAX) == 0) .* gen(vload, QMIN);
if any( abs( gen(vload, QG) - gen(vload, PG) .* Qlim ./ gen(vload, PMIN) ) > 1e-4 )
    errstr = sprintf('%s\n', ...
        'For a dispatchable load, PG and QG must be consistent', ...
        'with the power factor defined by PMIN and the Q limits.' );
    error(errstr);
end

% Find out which generators require additional linear constraints
% (as opposed to simple box constraints) on (Pg,Qg) to correctly
% model their PQ capability curves
ipqh = find( hasPQcap(gen, 'U') );
ipql = find( hasPQcap(gen, 'L') );

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
ipwl = find(gencost(:, MODEL) == PW_LINEAR);  %% piece-wise linear costs
nb = size(bus, 1);      %% number of buses
nl = size(branch, 1);   %% number of branches
ng = size(gen, 1);      %% number of dispatchable injections
nx = 2*nb + 2*ng;       %% number of standard OPF control variables
ny = size(ipwl, 1);     %% number of piece-wise linear costs
nvl  = size(vload, 1);  %% number of dispatchable loads
npqh = size(ipqh, 1);   %% number of general PQ capability curves (upper)
npql = size(ipql, 1);   %% number of general PQ capability curves (lower)
nusr = size(Au, 1);     %% number of linear user constraints
nw = size(N, 1);        %% number of general cost vars, w
if isempty(Au)
  nz = 0;
  Au = sparse(0,nx);
  if ~isempty(N)        %% still need to check number of columns of N
    if size(N, 2) ~= nx;
      error(sprintf('fmincopf.m: user supplied N matrix must have %d columns.', nx));
    end
  end
else
  nz = size(Au,2) - nx; %% additional linear variables
  if nz < 0
    error(sprintf('fmincopf.m: user supplied A matrix must have at least %d columns.', nx));
  end
end
nxyz = nx+ny+nz;        %% total number of control vars of all types

%% define indexing of optimization variable vector
k = 0;
thbas   = k + 1;    k = k + nb;     thend = k;      %% voltage angles
vbas    = k + 1;    k = k + nb;     vend  = k;      %% voltage magnitudes
pgbas   = k + 1;    k = k + ng;     pgend = k;      %% real power injections
qgbas   = k + 1;    k = k + ng;     qgend = k;      %% reactive power injections
ybas    = k + 1;    k = k + ny;     yend  = k;      %% pwl costs
zbas    = k + 1;    k = k + nz;     zend  = k;      %% user vars

%% define indexing of constraint vector
k = 0;
pmsmbas = k + 1;    k = k + nb;     pmsmend = k;    %% real power balance
qmsmbas = k + 1;    k = k + nb;     qmsmend = k;    %% reactive power balance
sfbas   = k + 1;    k = k + nl;     sfend   = k;    %% "from" flow limit
stbas   = k + 1;    k = k + nl;     stend   = k;    %% "to" flow limit
usrbas  = k + 1;    k = k + nusr;   usrend  = k;    % warning: nusr could be 0
pqhbas  = k + 1;    k = k + npqh;   pqhend  = k;    % warning: npqh could be 0
pqlbas  = k + 1;    k = k + npql;   pqlend  = k;    % warning: npql could be 0
vlbas   = k + 1;    k = k + nvl;    vlend   = k;    % warning: nvl could be 0
angbas  = k + 1;    k = k + nang;   angend  = k;    %% branch angle diff lims
%% not done yet, need number of Ay constraints.

% Let makeAy deal with any y-variable for piecewise-linear convex costs.
% note that if there are z variables then Ay doesn't have the columns
% that would span the z variables, so we append them.
if ny > 0
  [Ay, by]  = makeAy(baseMVA, ng, gencost, pgbas, qgbas, ybas);
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

% Make Avl, lvl, uvl in case there is a need for dispatchable loads
if nvl > 0
  xx = gen(vload, PMIN);
  yy = Qlim;
  pftheta = atan2(yy, xx);
  pc = sin(pftheta);
  qc = -cos(pftheta);
  ii = [ (1:nvl)'; (1:nvl)' ];
  jj = [ pgbas+vload-1; qgbas+vload-1 ];
  Avl = sparse(ii, jj, [pc; qc], nvl, nxyz);
  lvl = zeros(nvl, 1);
  uvl = lvl;
else
  Avl =[];
  lvl =[];
  uvl =[];
end

% Make Apqh if there is a need to add general PQ capability curves;
% use normalized coefficient rows so multipliers have right scaling
% in $$/pu
if npqh > 0
  Apqhdata = [gen(ipqh,QC1MAX)-gen(ipqh,QC2MAX), gen(ipqh,PC2)-gen(ipqh,PC1)];
  ubpqh = (gen(ipqh,QC1MAX)-gen(ipqh,QC2MAX)) .* gen(ipqh,PC1) ...
         + (gen(ipqh,PC2)-gen(ipqh,PC1)) .* gen(ipqh,QC1MAX);
  for i=1:npqh,
    tmp = norm(Apqhdata(i,:));
    Apqhdata(i,:) = Apqhdata(i, :) / tmp;
    ubpqh(i) = ubpqh(i) / tmp;
  end
  Apqh = sparse([1:npqh, 1:npqh]', [(pgbas-1)+ipqh;(qgbas-1)+ipqh], ...
                Apqhdata(:), npqh, nxyz);
  ubpqh = ubpqh / baseMVA;
  lbpqh = -1e10*ones(npqh,1);
else
  Apqh = [];
  ubpqh = [];
  lbpqh = [];
end

% similarly Apql
if npql > 0
  Apqldata = [gen(ipql,QC2MIN)-gen(ipql,QC1MIN), gen(ipql,PC1)-gen(ipql,PC2)];
  ubpql= (gen(ipql,QC2MIN)-gen(ipql,QC1MIN)) .* gen(ipql,PC1) ...
         - (gen(ipql,PC2)-gen(ipql,PC1)) .* gen(ipql,QC1MIN) ;
  for i=1:npql,
    tmp = norm(Apqldata(i, : ));
    Apqldata(i, :) = Apqldata(i, :) / tmp;
    ubpql(i) = ubpql(i) / tmp;
  end
  Apql = sparse([1:npql, 1:npql]', [(pgbas-1)+ipql;(qgbas-1)+ipql], ...
                Apqldata(:), npql, nxyz);
  ubpql = ubpql / baseMVA;
  lbpql = -1e10*ones(npql,1);
else
  Apql = [];
  ubpql = [];
  lbpql = [];
end

% % reorder columns of Au and N according to new gen ordering
% if ~isempty(Au)
%   if nz > 0
%     Au = Au(:, [[1:vend]'; vend+[igen; ng+igen]; qgend+[1:nz]']);
%   else
%     Au = Au(:, [[1:vend]'; vend+[igen; ng+igen]]);
%   end
% end
% if ~isempty(N)
%   if nz > 0
%     N =  N(:, [[1:vend]'; vend+[igen; ng+igen]; qgend+[1:nz]']);
%   else
%     N =  N(:, [[1:vend]'; vend+[igen; ng+igen]]);
%   end
% end

% Insert y columns in Au and N as necessary
if ny > 0
  if nz > 0
    Au = [ Au(:,1:qgend) sparse(nusr, ny) Au(:, qgend+(1:nz)) ];
    if ~isempty(N)
        N = [ N(:,1:qgend) sparse(nw, ny) N(:, qgend+(1:nz)) ];
    end
  else
    Au = [ Au sparse(nusr, ny) ];
    if ~isempty(N)
        N = [ N sparse(nw, ny) ];
    end
  end
end

% Now form the overall linear restriction matrix;
% note the order of the constraints.

if (ncony > 0 )
  A = [ Au;
        Apqh;
        Apql;
        Avl;
        Aang;
        Ay;
        sparse(ones(1,ny), ybas:yend, ones(1,ny), 1, nxyz ) ];  % "linear" cost
  l = [ lbu;
        lbpqh;
        lbpql;
        lvl;
        lang;
       -1e10*ones(ncony+1, 1) ];
  u = [ ubu;
        ubpqh;
        ubpql;
        uvl;
        uang;
        by;
        1e10];
else
  A = [ Au; Apqh; Apql; Avl; Aang ];
  l = [ lbu; lbpqh; lbpql; lvl; lang ];
  u = [ ubu; ubpqh; ubpql; uvl; uang ];
end


% So, can we do anything good about lambda initialization?
if all(bus(:, LAM_P) == 0)
  bus(:, LAM_P) = (10)*ones(nb, 1);
end


% --------------------------------------------------------------
% Up to this point, the setup is MINOS-like.  We now adapt
% things for fmincon.

% Form a vector with basic info to pass on as a parameter
parms = [ ...
    nb  ;% 1
    ng  ;% 2
    nl  ;% 3
    ny  ;% 4
    nx  ;% 5
    nvl ;% 6
    nz  ;% 7
    nxyz;% 8
    thbas;% 9
    thend;% 10
    vbas;% 11
    vend;% 12
    pgbas;% 13
    pgend;% 14
    qgbas;% 15
    qgend;% 16
    ybas;% 17
    yend;% 18
    zbas;% 19
    zend;% 20
    pmsmbas;% 21
    pmsmend;% 22
    qmsmbas;% 23
    qmsmend;% 24
    sfbas;% 25
    sfend;% 26
    stbas;% 27
    stend;% 28
];

% If there are y variables the last row of A is a linear cost vector
% of length nxyz. Let us excise it from A explicitly if it exists;
% otherwise it is zero.
if ny > 0
  nn = size(A,1);
  ccost = full(A(nn, :));
  A(nn, :) = [];
  l(nn) = [];
  u(nn) = [];
else
  ccost = zeros(1, nxyz);
end

% Divide l <= A*x <= u into less than, equal to, greater than, doubly-bounded
% sets.
ieq = find( abs(u-l) <= eps );          %% equality
igt = find( u >=  1e10 & l > -1e10 );   %% greater than, unbounded above
ilt = find( l <= -1e10 & u <  1e10 );   %% less than, unbounded below
ibx = find( (abs(u-l) > eps) & (u < 1e10) & (l > -1e10) );
Af  = [ A(ilt, :); -A(igt, :); A(ibx, :); -A(ibx, :) ];
bf  = [ u(ilt);   -l(igt);     u(ibx);    -l(ibx)];
Afeq = A(ieq, :);
bfeq = u(ieq);

% bounds on optimization vars; y vars unbounded
Va   = bus(:, VA) * (pi/180);
Pg   = gen(:, PG) / baseMVA;
Qg   = gen(:, QG) / baseMVA;
Pmin = gen(:, PMIN) / baseMVA;
Pmax = gen(:, PMAX) / baseMVA;
Qmin = gen(:, QMIN) / baseMVA;
Qmax = gen(:, QMAX) / baseMVA;
UB = Inf * ones(nxyz, 1);
LB = -UB;
LB(thbas-1+refs) = Va(refs);        UB(thbas-1+refs) = Va(refs);
LB(vbas:vend)    = bus(:, VMIN);    UB(vbas:vend)   = bus(:, VMAX);
LB(pgbas:pgend)  = Pmin;            UB(pgbas:pgend) = Pmax;
LB(qgbas:qgend)  = Qmin;            UB(qgbas:qgend) = Qmax;
if ~isempty(zl)
  LB(zbas:zend) = zl;
end
if ~isempty(zu)
  UB(zbas:zend) = zu;
end

% initialize optimization vars
x0 = zeros(nxyz, 1);
x0(thbas:thend) = Va;
x0(vbas:vend)   = bus(:, VM);
x0(vbas-1+gen(:,GEN_BUS)) = gen(:, VG);     % buses w/gens init V from gen data
x0(pgbas:pgend) = Pg;
x0(qgbas:qgend) = Qg;
% no ideas to initialize y variables
if ~isempty(z0)
  x0(zbas:zend) = z0;
end

% build admittance matrices
[Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);

% Tolerances
if mpopt(19) == 0   % # iterations
  mpopt(19) = 150 + 2*nb;
end

% basic optimset options needed for fmincon
fmoptions = optimset('GradObj', 'on', 'GradConstr', 'on' );
fmoptions = optimset(fmoptions, 'MaxIter', mpopt(19), 'TolCon', mpopt(16) );
fmoptions = optimset(fmoptions, 'TolX', mpopt(17), 'TolFun', mpopt(18) );
fmoptions.MaxFunEvals = 4 * fmoptions.MaxIter;
if verbose == 0,
  fmoptions.Display = 'off';
elseif verbose == 1
  fmoptions.Display = 'iter';
else
  fmoptions.Display = 'testing';
end

% select algorithm
otver = ver('optim');
if str2num(otver.Version(1)) < 4
  fmoptions = optimset(fmoptions, 'LargeScale', 'off');
  Af = full(Af);
  Afeq = full(Afeq);
else
  if mpopt(55) == 1           %% active-set
    fmoptions = optimset(fmoptions, 'Algorithm', 'active-set');
    Af = full(Af);
    Afeq = full(Afeq);
  elseif mpopt(55) == 2       %% interior-point, w/ default 'bfgs' Hessian approx
    fmoptions = optimset(fmoptions, 'Algorithm', 'interior-point');
  elseif mpopt(55) == 3       %% interior-point, w/ 'lbfgs' Hessian approx
    fmoptions = optimset(fmoptions, 'Algorithm', 'interior-point', 'Hessian','lbfgs');
  elseif mpopt(55) == 4       %% interior-point, w/ exact user-supplied Hessian
    fmc_hessian = @(x, lambda)hessfmin(x, lambda, baseMVA, bus, gen, gencost, branch, Ybus, Yf, Yt, mpopt, parms, ccost, N, fparm, H, Cw);
    fmoptions = optimset(fmoptions, 'Algorithm', 'interior-point', ...
        'Hessian', 'user-supplied', 'HessFcn', fmc_hessian);
  elseif mpopt(55) == 5       %% interior-point, w/ finite-diff Hessian
    fmoptions = optimset(fmoptions, 'Algorithm', 'interior-point', 'Hessian','fin-diff-grads', 'SubProblem', 'cg');
  else
    error('fmincopf: unknown algorithm specified in FMC_ALG option');
  end
end
% fmoptions = optimset(fmoptions, 'DerivativeCheck', 'on', 'FinDiffType', 'central', 'FunValCheck', 'on');
% fmoptions = optimset(fmoptions, 'Diagnostics', 'on');

if str2num(otver.Version(1)) >= 4 & strcmp(optimget(fmoptions, 'Algorithm'), 'interior-point')
  % set initial point
  x0 = zeros(nxyz, 1);
  x0(thbas:thend) = 0;
  x0(vbas:vend)   = 1;
  x0(pgbas:pgend) = (gen(:, PMIN) + gen(:, PMAX)) / 2 / baseMVA;
  x0(qgbas:qgend) = (gen(:, QMIN) + gen(:, QMAX)) / 2 / baseMVA;
  % no ideas to initialize y variables
  if ~isempty(zl) & ~isempty(zu)
    x0(zbas:zend) = (zl + zu) / 2;
  end
end

%%-----  run opf  -----
fmc_cost = @(x)costfmin(x, baseMVA, gencost, parms, ccost, N, fparm, H, Cw);
fmc_cons = @(x)consfmin(x, baseMVA, bus, gen, branch, Ybus, Yf, Yt, mpopt, parms);
[x, f, info, Output, Lambda, Jac] = ...
  fmincon(fmc_cost, x0, Af, bf, Afeq, bfeq, LB, UB, fmc_cons, fmoptions);

success = (info > 0);

% Unpack optimal x
Va = x(thbas:thend);
Vm = x(vbas:vend);
V = Vm .* exp(j*Va);
bus(:, VA) = Va * 180/pi;
bus(:, VM) = Vm;
gen(:, PG) = baseMVA * x(pgbas:pgend);
gen(:, QG) = baseMVA * x(qgbas:qgend);
gen(:, VG) = Vm(gen(:, GEN_BUS));

%% compute branch injections
Sf = V(branch(:, F_BUS)) .* conj(Yf * V);  %% cplx pwr at "from" bus, p.u.
St = V(branch(:, T_BUS)) .* conj(Yt * V);  %% cplx pwr at "to" bus, p.u.
branch(:, PF) = real(Sf) * baseMVA;
branch(:, QF) = imag(Sf) * baseMVA;
branch(:, PT) = real(St) * baseMVA;
branch(:, QT) = imag(St) * baseMVA;

%% line constraint is actually on square of limit
%% so we must fix multipliers
Lambda.ineqnonlin(1:nl)      = 2 * Lambda.ineqnonlin(1:nl)      .* branch(:, RATE_A) / baseMVA;
Lambda.ineqnonlin(nl+1:2*nl) = 2 * Lambda.ineqnonlin(nl+1:2*nl) .* branch(:, RATE_A) / baseMVA;

% Put in Lagrange multipliers
gen(:, MU_PMAX)  = Lambda.upper(pgbas:pgend) / baseMVA;
gen(:, MU_PMIN)  = Lambda.lower(pgbas:pgend) / baseMVA;
gen(:, MU_QMAX)  = Lambda.upper(qgbas:qgend) / baseMVA;
gen(:, MU_QMIN)  = Lambda.lower(qgbas:qgend) / baseMVA;
bus(:, LAM_P)    = Lambda.eqnonlin(1:nb) / baseMVA;
bus(:, LAM_Q)    = Lambda.eqnonlin(nb+1:2*nb) / baseMVA;
bus(:, MU_VMAX)  = Lambda.upper(vbas:vend);
bus(:, MU_VMIN)  = Lambda.lower(vbas:vend);
branch(:, MU_SF) = Lambda.ineqnonlin(1:nl) / baseMVA;
branch(:, MU_ST) = Lambda.ineqnonlin(nl+1:2*nl) / baseMVA;

% extract lambdas from linear constraints
nlt = length(ilt);
ngt = length(igt);
nbx = length(ibx);
lam = zeros(size(u));
lam(ieq) = Lambda.eqlin;
lam(ilt) = Lambda.ineqlin(1:nlt);
lam(igt) = -Lambda.ineqlin(nlt+[1:ngt]);
lam(ibx) = Lambda.ineqlin(nlt+ngt+[1:nbx]) - Lambda.ineqlin(nlt+ngt+nbx+[1:nbx]);

% stick in non-linear constraints too, so we can use the indexing variables
% we've defined, and negate so it looks like the pimul from MINOS
pimul = [
  -Lambda.eqnonlin(1:nb);
  -Lambda.eqnonlin(nb+1:2*nb);
  -Lambda.ineqnonlin(1:nl);
  -Lambda.ineqnonlin(nl+1:2*nl);
  -lam;
  -1;       %% dummy entry corresponding to linear cost row in A (in MINOS)
  Lambda.lower - Lambda.upper
];

% If we succeeded and there were generators with general PQ curve
% characteristics, this is the time to re-compute the multipliers,
% splitting any nonzero multiplier on one of the linear bounds among the
% Pmax, Pmin, Qmax or Qmin limits, producing one multiplier for a P limit and
% another for a Q limit. For upper Q limit, if we are neither at Pmin nor at 
% Pmax, the limit is taken as Pmin if the Qmax line's normal has a negative P
% component, Pmax if it has a positive P component. Messy but there really
% are many cases.  Remember multipliers in pimul() are negative.
muPmax = gen(:, MU_PMAX);
muPmin = gen(:, MU_PMIN);
if success & (npqh > 0)
%   gen(:, [MU_PMIN MU_PMAX MU_QMIN MU_QMAX])
  k = 1;
  for i = ipqh'
    if muPmax(i) > 0
      gen(i,MU_PMAX)=gen(i,MU_PMAX)-pimul(pqhbas+k-1)*Apqhdata(k,1)/baseMVA;
    elseif muPmin(i) > 0
      gen(i,MU_PMIN)=gen(i,MU_PMIN)+pimul(pqhbas+k-1)*Apqhdata(k,1)/baseMVA;
    else
      if Apqhdata(k, 1) >= 0
         gen(i,MU_PMAX)=gen(i,MU_PMAX)-pimul(pqhbas+k-1)*Apqhdata(k,1)/baseMVA;
      else
         gen(i,MU_PMIN)=gen(i,MU_PMIN)+pimul(pqhbas+k-1)*Apqhdata(k,1)/baseMVA;
      end
    end
    gen(i,MU_QMAX)=gen(i,MU_QMAX)-pimul(pqhbas+k-1)*Apqhdata(k,2)/baseMVA;
    k = k + 1;
  end
end

if success & (npql > 0)
%   gen(:, [MU_PMIN MU_PMAX MU_QMIN MU_QMAX])
  k = 1;
  for i = ipql'
    if muPmax(i) > 0
      gen(i,MU_PMAX)=gen(i,MU_PMAX)-pimul(pqlbas+k-1)*Apqldata(k,1)/baseMVA;
    elseif muPmin(i) > 0
      gen(i,MU_PMIN)=gen(i,MU_PMIN)+pimul(pqlbas+k-1)*Apqldata(k,1)/baseMVA;
    else
      if Apqldata(k,1) >= 0
        gen(i,MU_PMAX)=gen(i,MU_PMAX)-pimul(pqlbas+k-1)*Apqldata(k,1)/baseMVA;
      else
        gen(i,MU_PMIN)=gen(i,MU_PMIN)+pimul(pqlbas+k-1)*Apqldata(k,1)/baseMVA;
      end
    end
    gen(i,MU_QMIN)=gen(i,MU_QMIN)+pimul(pqlbas+k-1)*Apqldata(k,2)/baseMVA;
    k = k + 1;
  end
%   gen(:, [MU_PMIN MU_PMAX MU_QMIN MU_QMAX])
%   -[ pimul(pqlbas-1+[1:2]) pimul(pqhbas-1+[1:2]) ]/baseMVA
%   -[ pimul(pqlbas-1+[1:2]).*Apqldata([1:2],1) pimul(pqhbas-1+[1:2]).*Apqhdata([1:2],1) ]/baseMVA
%   -[ pimul(pqlbas-1+[1:2]).*Apqldata([1:2],2) pimul(pqhbas-1+[1:2]).*Apqhdata([1:2],2) ]/baseMVA
end

% angle limit constraints
if success & (nang > 0)
  tmp = [angbas:angend];
  ii = find(pimul(tmp) > 0);
  branch(iang(ii), MU_ANGMIN) = pimul(tmp(ii)) * pi/180;
  ii = find(pimul(tmp) < 0);
  branch(iang(ii), MU_ANGMAX) = -pimul(tmp(ii)) * pi/180;
end

% With these modifications, printpf must then look for multipliers
% if available in order to determine if a limit has been hit.

% We are done with standard opf but we may need to provide the
% constraints and their Jacobian also.
if nargout > 7
  [g, geq, dg, dgeq] = consfmin(x, baseMVA, bus, gen, branch, ...
                                Ybus, Yf, Yt, mpopt, parms);
  g = [ geq; g];
  jac = [ dgeq'; dg']; % true Jacobian organization
end

% Go back to original data.
% reorder generators
gen = gen(inv_gen_ord, :);

% convert to original external bus ordering
[bus, gen, branch, areas] = int2ext(i2e, bus, gen, branch, areas);

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
if (nargout == 0) & ( success )
  printpf(baseMVA, bus, genout, branchout, f, info, et, 1, mpopt);
end

if nargout, busout = bus; end

return;

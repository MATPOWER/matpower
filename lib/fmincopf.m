function [busout, genout, branchout, f, success, info, et, g, jac] = ...
      fmincopf(baseMVA, bus, gen, branch, areas, gencost, Au, lbu, ubu, mpopt, ...
           N, fparm, H, Cw)
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
%   [bus, gen, branch, f, success, info, et, g, jac] = fmincopf(casefile)
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
%   fparm matrix (see manual or type 'help generalcost').  If the
%   resulting vector is now named w, then H and Cw define a quadratic
%   cost on w:  (1/2)*w'*H*w + Cw * w . H and N should be sparse matrices
%   and H should also be symmetric.
%
%   The optional mpopt vector specifies MATPOWER options. Type 'help mpoption'
%   for details and default values.
%
%   The solved case is returned in the data matrices, bus, gen and branch. Also
%   returned are the final objective function value (f) and a flag which is
%   true if the algorithm was successful in finding a solution (success).
%   Additional optional return values are an algorithm specific return status
%   (info), elapsed time in seconds (et), the constraint vector (g) and the
%   Jacobian matrix (jac).
%
%   Rules for A matrix: If the user specifies an A matrix that has more columns
%   than the number of "x" (OPF) variables, then there are extra linearly
%   constrained "z" variables.
%
%   NOTE: The shadow prices (lambda's and mu's) produced by fmincon appear to
%         be slightly inaccurate.

%   MATPOWER
%   $Id$
%   by Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Autonoma de Manizales
%   Copyright (c) 2000-2006 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

% Sort out input arguments
t1 = clock;
if isstr(baseMVA) | isstruct(baseMVA)   % passing filename or struct
  % 14  fmincopf(baseMVA,  bus, gen, branch, areas, gencost, Au,    lbu, ubu, mpopt, N, fparm, H, Cw)
  % 9   fmincopf(casefile, Au,  lbu, ubu,    mpopt, N,       fparm, H,   Cw)
  % 5   fmincopf(casefile, Au,  lbu, ubu,    mpopt)
  % 4   fmincopf(casefile, Au,  lbu, ubu)
  % 2   fmincopf(casefile, mpopt)
  % 1   fmincopf(casefile)
  if any(nargin == [1, 2, 4, 5, 9])
    casefile = baseMVA;
    if nargin == 9
      N     = gencost;
      fparm = Au;
      H     = lbu;
      Cw    = ubu;
    else
      N     = [];
      fparm = [];
      H     = [];
      Cw    = [];
    end
    if nargin < 4
      Au  = sparse(0,0);
      lbu = [];
      ubu = [];
    else
      Au  = bus;
      lbu = gen;
      ubu = branch;
    end
    if nargin == 9 | nargin == 5
      mpopt = areas;
    elseif nargin == 2
      mpopt = bus;
    else
      mpopt = [];
    end
  else
    error('fmincopf.m: Incorrect input parameter order, number or type');
  end
  [baseMVA, bus, gen, branch, areas, gencost] = loadcase(casefile);
else    % passing individual data matrices
  % 14  fmincopf(baseMVA,  bus, gen, branch, areas, gencost, Au,    lbu, ubu, mpopt, N, fparm, H, Cw)
  % 10  fmincopf(baseMVA,  bus, gen, branch, areas, gencost, Au,    lbu, ubu, mpopt)
  % 9   fmincopf(baseMVA,  bus, gen, branch, areas, gencost, Au,    lbu, ubu)
  % 7   fmincopf(baseMVA,  bus, gen, branch, areas, gencost, mpopt)
  % 6   fmincopf(baseMVA,  bus, gen, branch, areas, gencost)
  if any(nargin == [6, 7, 9, 10, 14])
    if nargin ~= 14
      N     = [];
      fparm = [];
      H     = [];
      Cw    = [];
    end
    if nargin == 7
      mpopt = Au;
    elseif nargin == 6 | nargin == 9
      mpopt = [];
    end
    if nargin < 9
      Au  = sparse(0,0);
      lbu = [];
      ubu = [];
    end
  else
    error('fmincopf.m: Incorrect input parameter order, number or type');
  end
end
if size(N, 1) > 0
  if size(N, 1) ~= size(fparm, 1) | size(N, 1) ~= size(H, 1) | ...
     size(N, 1) ~= size(H, 2) | size(N, 1) ~= length(Cw)
    error('fmincopf.m: wrong dimensions in generalized cost parameters');
  end
  if size(Au, 1) > 0 & size(N, 2) ~= size(Au, 2)
    error('fmincopf.m: A and N must have the same number of columns');
  end
end
if isempty(mpopt)
  mpopt = mpoption;
end

% Load column indexes for case tables.
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
[ref, pv, pq] = bustypes(bus, gen);

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
if size(find(bus(:, BUS_TYPE) == REF), 1) > 1
  errstr = ['\nfmincopf: Warning: more than one reference bus detected in bus table data.\n', ...
              '      For a system with islands, a reference bus in each island\n', ...
              '      might help convergence but in a fully connected system such\n', ...
              '      a situation is probably not reasonable.\n\n' ];
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


% Find out problem dimensions
nb = size(bus, 1);                              % buses
ng = size(gen, 1);                              % variable injections
nl = size(branch, 1);                           % branches
iycost = find(gencost(:, MODEL) == PW_LINEAR);  % y variables for pwl cost
ny    = size(iycost, 1);
neqc  = 2 * nb;                                 % nonlinear equalities
nusr  = size(Au, 1);                            % # linear user constraints
nx    = 2*nb + 2*ng;                            % control variables
nvl   = size(vload, 1);                         % dispatchable loads
npqh  = size(ipqh, 1);                          % general pq capability curves
npql  = size(ipql, 1);
if isempty(Au)
  nz = 0;
  Au = sparse(0,nx);
  if ~isempty(N)        % still need to check number of columns of N
    if size(N,2) ~= nx;
      error(sprintf('fmincopf.m: user supplied N matrix must have %d columns.', nx));
    end
  end
else
  nz = size(Au,2) - nx;                       % additional linear variables
  if nz < 0
    error(sprintf('fmincopf.m: user supplied A matrix must have at least %d columns.', nx));
  end
end
nxyz = nx+ny+nz;                                % total # of vars of all types

% Definition of indexes into optimization variable vector and constraint
% vector.
thbas = 1;                thend    = thbas+nb-1;
vbas     = thend+1;       vend     = vbas+nb-1;
pgbas    = vend+1;        pgend    = pgbas+ng-1;
qgbas    = pgend+1;       qgend    = qgbas+ng-1;
ybas     = qgend + 1;     yend     = ybas + ny - 1;
zbas     = yend + 1;      zend     = zbas + nz - 1;

pmsmbas = 1;              pmsmend = pmsmbas+nb-1;
qmsmbas = pmsmend+1;      qmsmend = qmsmbas+nb-1;
sfbas   = qmsmend+1;      sfend   = sfbas+nl-1;
stbas   = sfend+1;        stend   = stbas+nl-1;
usrbas  = stend+1;        usrend  = usrbas+nusr-1; % warning: nusr could be 0
pqhbas  = usrend+1;       pqhend  = pqhbas+npqh-1; % warning: npqh could be 0
pqlbas  = pqhend+1;       pqlend  = pqlbas+npql-1; % warning: npql could be 0
vlbas   = pqlend+1;       vlend   = vlbas+nvl-1; % not done yet, need # of
                                                 % Ay constraints.

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
yconbas = vlend+1;        yconend = yconbas+ncony-1; % finally done with
                                                     % constraint indexing

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

% Insert y columns in Au and N as necessary
if ny > 0
  if nz > 0
    Au = [ Au(:,1:qgend) sparse(nusr, ny) Au(:, qgend+(1:nz)) ];
    if ~isempty(N)
        N = [ N(:,1:qgend) sparse(size(N,1), ny) N(:, qgend+(1:nz)) ];
    end
  else
    Au = [ Au sparse(nusr, ny) ];
    if ~isempty(N)
        N = [ N sparse(size(N,1), ny) ];
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
        Ay;
        sparse(ones(1,ny), ybas:yend, ones(1,ny), 1, nxyz ) ];  % "linear" cost
  l = [ lbu;
        lbpqh;
        lbpql;
        lvl;
       -1e10*ones(ncony+1, 1) ];
  u = [ ubu;
        ubpqh;
        ubpql;
        uvl;
        by;
        1e10];
else
  A = [ Au; Apqh; Apql; Avl ];
  l = [ lbu; lbpqh; lbpql; lvl ];
  u = [ ubu; ubpqh; ubpql; uvl ];
end


% So, can we do anything good about lambda initialization?
if all(bus(:, LAM_P) == 0)
  bus(:, LAM_P) = (10)*ones(nb, 1);
end


% --------------------------------------------------------------
% Up to this point, the setup is MINOS-like.  We now adapt
% things for fmincon.

j = sqrt(-1);

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
ieq = find( abs(u-l) <= eps);
igt = find( u >= 1e10);  % unlimited ceiling
ilt = find( l <= -1e10); % unlimited bottom
ibx = find( (abs(u-l) > eps) & (u < 1e10) & (l > -1e10));
Af  = [ A(ilt, :); -A(igt, :); A(ibx, :); -A(ibx, :) ];
bf  = [ u(ilt);   -l(igt);     u(ibx);    -l(ibx)];
Afeq = A(ieq, :);
bfeq = u(ieq);

% bounds on optimization vars; y and z vars unbounded at box bounds;
% if needed, user must do this via the Ax < l mechanism if needed and hope that
% fmincon handles singleton rows elegantly.
UB = Inf * ones(nxyz, 1);
LB = -UB;
LB(thbas+ref-1) = bus(ref, VA)*pi/180;  UB(thbas+ref-1) = bus(ref, VA)*pi/180;
LB(vbas:vend)   = bus(:, VMIN);         UB(vbas:vend)   = bus(:, VMAX);
LB(pgbas:pgend) = gen(:, PMIN)/baseMVA; UB(pgbas:pgend) = gen(:, PMAX)/baseMVA;
LB(qgbas:qgend) = gen(:, QMIN)/baseMVA; UB(qgbas:qgend) = gen(:, QMAX)/baseMVA;

% Compute initial vector
x0 = zeros(nxyz, 1);
x0(thbas:thend) = bus(:, VA) * pi/180;
x0(vbas:vend)   = bus(:, VM);
x0(vbas+gen(:,GEN_BUS)-1) = gen(:, VG);   % buses w. gens init V from gen data
x0(pgbas:pgend) = gen(:, PG) / baseMVA;
x0(qgbas:qgend) = gen(:, QG) / baseMVA;
% no ideas to initialize z, y variables, though, and no mechanism yet
% to ask for user-provided initial z, y.

% build admittance matrices
[Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);

% Tolerances
if mpopt(19) == 0   % # iterations
  mpopt(19) = 150 + 2*nb;
end

% basic optimset options needed for fmincon
% fmoptions = optimset('GradObj', 'on', 'Hessian', 'on', 'LargeScale', 'on', ...
%                    'GradConstr', 'on');
fmoptions = optimset('GradObj', 'on', 'LargeScale', 'off', 'GradConstr', 'on');
fmoptions = optimset(fmoptions, 'MaxIter', mpopt(19), 'TolCon', mpopt(16) );
fmoptions = optimset(fmoptions, 'TolX', mpopt(17), 'TolFun', mpopt(18) );
fmoptions.MaxFunEvals = 4 * fmoptions.MaxIter;
if mpopt(31) == 0,
  fmoptions.Display = 'off';
else
  fmoptions.Display = 'iter';
end

Af = full(Af);
Afeq = full(Afeq);
[x, f, info, Output, Lambda, Jac] = ...
  fmincon('costfmin', x0, Af, bf, Afeq, bfeq, LB, UB, 'consfmin', fmoptions, ...
         baseMVA, bus, gen, gencost, branch, areas, Ybus, Yf, Yt, mpopt, ...
         parms, ccost, N, fparm, H, Cw);
success = (info > 0);

% Unpack optimal x
bus(:, VA) = x(thbas:thend)*180/pi;
bus(:, VM) = x(vbas:vend);
gen(:, PG) = baseMVA * x(pgbas:pgend);
gen(:, QG) = baseMVA * x(qgbas:qgend);
gen(:, VG) = bus(gen(:, GEN_BUS), VM);

% reconstruct voltages
Va = x(thbas:thend);
Vm = x(vbas:vend);
V = Vm .* exp(j*Va);

%% compute branch injections
Sf = V(branch(:, F_BUS)) .* conj(Yf * V);  %% cplx pwr at "from" bus, p.u.
St = V(branch(:, T_BUS)) .* conj(Yt * V);  %% cplx pwr at "to" bus, p.u.
branch(:, PF) = real(Sf) * baseMVA;
branch(:, QF) = imag(Sf) * baseMVA;
branch(:, PT) = real(St) * baseMVA;
branch(:, QT) = imag(St) * baseMVA;

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
lam(igt) = Lambda.ineqlin(nlt+[1:ngt]);
lam(ibx) = Lambda.ineqlin(nlt+ngt+[1:nbx]) + Lambda.ineqlin(nlt+ngt+nbx+[1:nbx]);

% stick in non-linear constraints too, so we can use the indexing variables
% we've defined, and negate so it looks like the pimul from MINOS
pimul = [ zeros(stend,1); -lam ];

% If we succeeded and there were generators with general pq curve
% characteristics, this is the time to re-compute the multipliers,
% splitting any nonzero multiplier on one of the linear bounds among the
% Pmax, Pmin, Qmax or Qmin limits, producing one multiplier for a P limit and
% another for a Q limit. For upper Q limit, if we are neither at Pmin nor at 
% Pmax, the limit is taken at Pmin if the Qmax line's normal has a negative P
% component, Pmax if it has a positive P component. Messy but there really
% are many cases.
if success & (npqh > 0)
  k = 1;
  for i = ipqh'
    if gen(i, MU_PMAX) > 0
      gen(i,MU_PMAX)=gen(i,MU_PMAX)-pimul(pqhbas+k-1)*Apqhdata(k,1)/baseMVA;
      gen(i,MU_QMAX)=gen(i,MU_QMAX)-pimul(pqhbas+k-1)*Apqhdata(k,2)/baseMVA;
    elseif gen(i, MU_PMIN) < 0
      gen(i,MU_PMIN)=gen(i,MU_PMIN)+pimul(pqhbas+k-1)*Apqhdata(k,1)/baseMVA;
      gen(i,MU_QMAX)=gen(i,MU_QMAX)-pimul(pqhbas+k-1)*Apqhdata(k,2)/baseMVA;
    else
      if Apqhdata(k, 1) >= 0
         gen(i, MU_PMAX) = -pimul(pqhbas+k-1)*Apqhdata(k,1)/baseMVA;
      else
         gen(i, MU_PMIN) = pimul(pqhbas+k-1)*Apqhdata(k,1)/baseMVA;
      end
      gen(i, MU_QMAX)= gen(i,MU_QMAX)-pimul(pqhbas+k-1)*Apqhdata(k,2)/baseMVA;
    end
    k = k + 1;
  end
end

if success & (npql > 0)
  k = 1;
  for i = ipql'
    if gen(i, MU_PMAX) > 0
      gen(i,MU_PMAX)=gen(i,MU_PMAX)-pimul(pqlbas+k-1)*Apqldata(k,1)/baseMVA;
      gen(i,MU_QMIN)=gen(i,MU_QMIN)-pimul(pqlbas+k-1)*Apqldata(k,2)/baseMVA;
    elseif gen(i, MU_PMIN) > 0
      gen(i,MU_PMIN)=gen(i,MU_PMIN)-pimul(pqlbas+k-1)*Apqldata(k,1)/baseMVA;
      gen(i,MU_QMIN)=gen(i,MU_QMIN)+pimul(pqlbas+k-1)*Apqldata(k,2)/baseMVA;
    else
      if Apqldata(k,1) >= 0
        gen(i,MU_PMAX)= -pimul(pqlbas+k-1)*Apqldata(k,1)/baseMVA;
      else
        gen(i,MU_PMIN)= pimul(pqlbas+k-1)*Apqldata(k,1)/baseMVA;
      end
      gen(i,MU_QMIN)=gen(i,MU_QMIN)+pimul(pqlbas+k-1)*Apqldata(k,2)/baseMVA;
    end
    k = k + 1;
  end
end

% With these modifications, printpf must then look for multipliers
% if available in order to determine if a limit has been hit.

% We are done with standard opf but we may need to provide the
% constraints and their Jacobian also.
if nargout > 7
  [g, geq, dg, dgeq] = consfmin(x, baseMVA, bus, gen, gencost, branch, areas,...
                                Ybus, Yf, Yt, mpopt, parms, ccost);
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

et = etime(clock,t1);
if  (nargout == 0) & ( success )
  printpf(baseMVA, bus, genout, branchout, f, info, et, 1, mpopt);
end

if nargout, busout = bus; end

return;

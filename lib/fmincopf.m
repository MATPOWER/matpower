function [busout, genout, branchout, f, success, info, et, g, jac] = ...
      fmincopf(baseMVA, bus, gen, branch, areas, gencost, Au, lbu, ubu, mpopt)
%FMINCOPF  Solves an AC optimal power flow using FMINCON (Opt Tbx 2.x & later).
%
%   [bus, gen, branch, f, success] = fmincopf(casefile, mpopt)
%
%   [bus, gen, branch, f, success] = fmincopf(casefile, A, l, u, mpopt)
%
%   [bus, gen, branch, f, success] = fmincopf(baseMVA, bus, gen, branch,...
%                                    areas,  gencost, mpopt)
%
%   [bus, gen, branch, f, success] = fmincopf(baseMVA, bus, gen, branch,...
%                                    areas, gencost, A, l, u, mpopt)
%
%   [bus, gen, branch, f, success, info, et, g, jac] = fmincopf(casefile)
%
%   The data for the problem can be specified in one of 3 ways: (1) the name of
%   a case file which defines the data matrices baseMVA, bus, gen, branch,
%   areas and gencost, (2) a struct containing the data matrices as fields, or
%   (3) the data matrices themselves.
%
%   When specified, A, l, u represent additional linear constraints on the
%   optimization variables, l <= A*x <= u. For an explanation of the
%   formulation used and instructions for forming the A matrix, type
%   'help genform'.
%
%   The optional mpopt vector specifies MATPOWER options. Type 'help mpoption'
%   for details and default values.
%
%   The solved case is returned in the data matrices, bus, gen and branch. Also,
%   returned are the final objective function value (f) and a flag which is
%   true if the algorithm was successful in finding a solution (success).
%   Additional optional return values are an algorithm specific return status
%   (info), elapsed time in seconds (et), the constraint vector (g) and the
%   Jacobian matrix (jac).
% 
%   Rules for A matrix: If the user specifies an A matrix that has more columns
%   than the combined number of "x" (OPF) and "y" (pwl cost) variables, thus
%   allowing for extra linearly costed and linearly constrained "z" variables,
%   then it is the user's responsibility to fill the overall linear cost in the
%   last row of A, including that reflected in the "y" variables.
%
%   NOTE: The shadow prices (lambda's and mu's) produced by fmincon appear to
%         be slightly inaccurate.

%   MATPOWER
%   $Id$
%   by Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Autonoma de Manizales
%   Copyright (c) 2000-2005 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

% Sort out input arguments
t1 = clock;
if isstr(baseMVA) | isstruct(baseMVA)
  casefile = baseMVA;
  if nargin == 5
    Au = bus;
    lbu = gen;
    ubu = branch;
    mpopt = areas;
  elseif nargin == 4
    Au = bus;
    lbu = gen;
    ubu = branch;
    mpopt = mpoption;
  elseif nargin == 2
    Au = sparse(0,0);
    lbu = [];
    ubu = [];
    mpopt = bus;
  elseif nargin == 1
    Au = sparse(0,0);
    lbu = [];
    ubu = [];
    mpopt = mpoption;
  else
    error('fmincopf.m: Incorrect input parameter order, number or type');
  end;
  [baseMVA, bus, gen, branch, areas, gencost] = loadcase(casefile);
else
  if nargin == 9 
    mpopt = mpoption;
  elseif nargin == 7
    mpopt = Au;
    Au = sparse(0,0);
    lbu = [];
    ubu = [];
  elseif nargin == 6
    mpopt = mpoption;
    Au = sparse(0,0);
    lbu = [];
    ubu = [];
  elseif nargin ~= 10
    error('fmincopf.m: Incorrect input parameter order, number or type');
  end
end

% Load column indexes for case tables.
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, ...
    PMAX, PMIN, MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, QMAX2, QMIN2, ...
    RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;

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

% Find out if any of these "generators" are actually dispatchable loads.
% (see 'help isload' for details on what consitutes a dispatchable load)
% Dispatchable loads are modeled as generators with an added constant
% power factor constraint. The power factor is derived from the original
% value of Pmin and either Qmin (for inductive loads) or Qmax (for capacitive
% loads). If both Qmin and Qmax are zero, this implies a unity power factor
% without the need for an additional constraint.
vload = find( isload(gen) & (gen(:, QMIN) ~= 0 | gen(:, QMAX) ~= 0) );
% At least one of the Q limits must be zero (corresponding to Pmax == 0)
if any( gen(vload, QMIN) ~= 0 & gen(vload, QMAX) ~= 0 )
	error('Either Qmin or Qmax must be equal to zero for each dispatchable load.');
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

% Find out problem dimensions
nb = size(bus, 1);                             % buses
ng = size(gen, 1);                             % variable injections
nl = size(branch, 1);                          % branches
iycost = find(gencost(:, MODEL) == PW_LINEAR); % y variables for pwl cost
ny = size(iycost, 1);
neqc = 2 * nb;                                 % nonlinear equalities
nx = 2*nb + 2*ng;                              % control variables
nvl = size(vload, 1);                          % dispatchable loads
nz = size(Au,2) - 2*nb - 2*ng - ny;            % number of extra z variables
nz = max(nz,0);

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

% Let makeAy deal with any y-variable for piecewise-linear convex costs.
[Ay, by]  = makeAy(baseMVA, ng, gencost, pgbas, qgbas, ybas);
ncony = size(Ay,1);

% Make Avl, lvl, uvl in case there is a need for dispatchable loads
if nvl > 0
  xx = gen(vload, PMIN);
  yy = Qlim;
  pftheta = atan2(yy, xx);
  pc = sin(pftheta);
  qc = -cos(pftheta);
  ii = [ (1:nvl)'; (1:nvl)' ];
  jj = [ pgbas+vload-1; qgbas+vload-1 ];
  Avl = sparse(ii, jj, [pc; qc], nvl, yend);
  lvl = zeros(nvl, 1);
  uvl = lvl;
else
  Avl =[];
  lvl =[];
  uvl =[];
end

% Now form the overall linear restriction matrix; note the order
% of the constraints.

if (nz > 0)
  % user defined z variables thus becoming responsible for
  % defining any Ay needed as well as the cost row (whether there are
  % y variables or not) in the last row of Au.
  A = [ Avl  sparse(size(Avl,1), size(Au,2)-size(Avl,2));
        Au ];
  l = [ lvl;
        lbu ];
  u = [ uvl;
        ubu ];
else               % no z variables
  if (ncony > 0 )  % ... but some y variables from pwl costs; we supply linear
    A = [ Au;      % cost for y variables in last row.
          Avl;
          Ay; 
        sparse(ones(1,ny), ybas:yend, ones(1,ny), 1, yend ) ]; % "linear" cost
    l = [ lbu; 
          lvl;
         -1e10*ones(ncony+1, 1) ];
    u = [ ubu;
          uvl;
          by;
          1e10];
  else                 % No y variables (no pwl costs) and no z variables
    A = [ Au;  Avl ];  % but perhaps we have user linear constraints in Au
    l = [ lbu; lvl ];  % on (theta,V,Pg,Qg) variables.
    u = [ ubu; uvl ];
  end
end

% So, can we do anything good about lambda initialization?
if all(bus(:, LAM_P) == 0)
  bus(:, LAM_P) = (10)*ones(nb, 1);
end

% total number of variables
nxyz = nx+ny+nz;


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

% If there are z variables or y variables the last row of A should be
% holding a linear cost vector of length nxyz.  Let us excise it from A
% explicitly if it exists; otherwise it is zero.
if ny+nz > 0
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
         parms, ccost);
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

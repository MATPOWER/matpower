function [busout, genout, branchout, f, success, info, et, g, jac, xr, pimul] = ...
    opf(varargin)
%OPF  Solves an optimal power flow.
%   [RESULTS, SUCCESS] = OPF(MPC, MPOPT)
%
%   Returns either a RESULTS struct and an optional SUCCESS flag, or individual
%   data matrices, the objective function value and a SUCCESS flag. In the
%   latter case, there are additional optional return values. See Examples
%   below for the possible calling syntax options.
%
%   Examples:
%       Output argument options:
%
%       results = opf(...)
%       [results, success] = opf(...)
%       [bus, gen, branch, f, success] = opf(...)
%       [bus, gen, branch, f, success, info, et, g, jac, xr, pimul] = opf(...)
%
%       Input arguments options:
%
%       opf(mpc)
%       opf(mpc, mpopt)
%       opf(mpc, userfcn, mpopt)
%       opf(mpc, A, l, u)
%       opf(mpc, A, l, u, mpopt)
%       opf(mpc, A, l, u, mpopt, N, fparm, H, Cw)
%       opf(mpc, A, l, u, mpopt, N, fparm, H, Cw, z0, zl, zu)
%
%       opf(baseMVA, bus, gen, branch, areas, gencost)
%       opf(baseMVA, bus, gen, branch, areas, gencost, mpopt)
%       opf(baseMVA, bus, gen, branch, areas, gencost, userfcn, mpopt)
%       opf(baseMVA, bus, gen, branch, areas, gencost, A, l, u)
%       opf(baseMVA, bus, gen, branch, areas, gencost, A, l, u, mpopt)
%       opf(baseMVA, bus, gen, branch, areas, gencost, A, l, u, ...
%                                   mpopt, N, fparm, H, Cw)
%       opf(baseMVA, bus, gen, branch, areas, gencost, A, l, u, ...
%                                   mpopt, N, fparm, H, Cw, z0, zl, zu)
%
%   The data for the problem can be specified in one of three ways:
%   (1) a string (mpc) containing the file name of a MATPOWER case
%     which defines the data matrices baseMVA, bus, gen, branch, and
%     gencost (areas is not used at all, it is only included for
%     backward compatibility of the API).
%   (2) a struct (mpc) containing the data matrices as fields.
%   (3) the individual data matrices themselves.
%   
%   The optional user parameters for user constraints (A, l, u), user costs
%   (N, fparm, H, Cw), user variable initializer (z0), and user variable
%   limits (zl, zu) can also be specified as fields in a case struct,
%   either passed in directly or defined in a case file referenced by name.
%   
%   When specified, A, l, u represent additional linear constraints on the
%   optimization variables, l <= A*[x; z] <= u. If the user specifies an A
%   matrix that has more columns than the number of "x" (OPF) variables,
%   then there are extra linearly constrained "z" variables. For an
%   explanation of the formulation used and instructions for forming the
%   A matrix, see the manual.
%
%   A generalized cost on all variables can be applied if input arguments
%   N, fparm, H and Cw are specified.  First, a linear transformation
%   of the optimization variables is defined by means of r = N * [x; z].
%   Then, to each element of r a function is applied as encoded in the
%   fparm matrix (see manual). If the resulting vector is named w,
%   then H and Cw define a quadratic cost on w: (1/2)*w'*H*w + Cw * w .
%   H and N should be sparse matrices and H should also be symmetric.
%
%   The optional mpopt vector specifies MATPOWER options. If the OPF
%   algorithm is not explicitly set in the options MATPOWER will use
%   the default solver, based on a primal-dual interior point method.
%   For the AC OPF this is OPF_ALG = 560, unless the TSPOPF optional
%   package is installed, in which case the default is 540. For the
%   DC OPF, the default is OPF_ALG_DC = 200. See MPOPTION for
%   more details on the available OPF solvers and other OPF options
%   and their default values.
%
%   The solved case is returned either in a single results struct (described
%   below) or in the individual data matrices, bus, gen and branch. Also
%   returned are the final objective function value (f) and a flag which is
%   true if the algorithm was successful in finding a solution (success).
%   Additional optional return values are an algorithm specific return status
%   (info), elapsed time in seconds (et), the constraint vector (g), the
%   Jacobian matrix (jac), and the vector of variables (xr) as well 
%   as the constraint multipliers (pimul).
%
%   The single results struct is a MATPOWER case struct (mpc) with the
%   usual baseMVA, bus, branch, gen, gencost fields, along with the
%   following additional fields:
%
%       .order      see 'help ext2int' for details of this field
%       .et         elapsed time in seconds for solving OPF
%       .success    1 if solver converged successfully, 0 otherwise
%       .om         OPF model object, see 'help opf_model'
%       .x          final value of optimization variables (internal order)
%       .f          final objective function value
%       .mu         shadow prices on ...
%           .var
%               .l  lower bounds on variables
%               .u  upper bounds on variables
%           .nln
%               .l  lower bounds on non-linear constraints
%               .u  upper bounds on non-linear constraints
%           .lin
%               .l  lower bounds on linear constraints
%               .u  upper bounds on linear constraints
%       .g          (optional) constraint values
%       .dg         (optional) constraint 1st derivatives
%       .df         (optional) obj fun 1st derivatives (not yet implemented)
%       .d2f        (optional) obj fun 2nd derivatives (not yet implemented)
%       .raw        raw solver output in form returned by MINOS, and more
%           .xr     final value of optimization variables
%           .pimul  constraint multipliers
%           .info   solver specific termination code
%           .output solver specific output information
%              .alg algorithm code of solver used
%       .var        
%           .val    optimization variable values, by named block
%               .Va     voltage angles
%               .Vm     voltage magnitudes (AC only)
%               .Pg     real power injections
%               .Qg     reactive power injections (AC only)
%               .y      constrained cost variable (only if have pwl costs)
%               (other) any user defined variable blocks
%           .mu     variable bound shadow prices, by named block
%               .l  lower bound shadow prices
%                   .Va, Vm, Pg, Qg, y, (other)
%               .u  upper bound shadow prices
%                   .Va, Vm, Pg, Qg, y, (other)
%       .nln    (AC only)
%           .mu     shadow prices on non-linear constraints, by named block
%               .l  lower bounds
%                   .Pmis   real power mismatch equations
%                   .Qmis   reactive power mismatch equations
%                   .Sf     flow limits at "from" end of branches
%                   .St     flow limits at "to" end of branches
%               .u  upper bounds
%                   .Pmis, Qmis, Sf, St
%       .lin
%           .mu     shadow prices on linear constraints, by named block
%               .l  lower bounds
%                   .Pmis   real power mistmatch equations (DC only)
%                   .Pf     flow limits at "from" end of branches (DC only)
%                   .Pt     flow limits at "to" end of branches (DC only)
%                   .PQh    upper portion of gen PQ-capability curve (AC only)
%                   .PQl    lower portion of gen PQ-capability curve (AC only)
%                   .vl     constant power factor constraint for loads (AC only)
%                   .ycon   basin constraints for CCV for pwl costs
%                   (other) any user defined constraint blocks
%               .u  upper bounds
%                   .Pmis, Pf, Pf, PQh, PQl, vl, ycon, (other)
%       .cost       user defined cost values, by named block
%
%   See also RUNOPF, DCOPF, UOPF, CASEFORMAT.

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   and Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Autonoma de Manizales
%   Copyright (c) 1996-2010 by Power System Engineering Research Center (PSERC)
%
%   This file is part of MATPOWER.
%   See http://www.pserc.cornell.edu/matpower/ for more info.
%
%   MATPOWER is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published
%   by the Free Software Foundation, either version 3 of the License,
%   or (at your option) any later version.
%
%   MATPOWER is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with MATPOWER. If not, see <http://www.gnu.org/licenses/>.
%
%   Additional permission under GNU GPL version 3 section 7
%
%   If you modify MATPOWER, or any covered work, to interface with
%   other modules (such as M-files and MEX-files) available in a
%   Matlab (or compatible) environment containing parts covered
%   under other licensing terms, the licensors of MATPOWER grant
%   you additional permission to convey the resulting work.

%%----- initialization -----
t0 = clock;         %% start timer

%% process input arguments
[mpc, mpopt] = opf_args(varargin{:});

%% options
dc  = mpopt(10);        %% PF_DC        : 1 = DC OPF, 0 = AC OPF 
alg = mpopt(11);        %% OPF_ALG
verbose = mpopt(31);    %% VERBOSE

%% set AC OPF algorithm code
if ~dc
  if alg == 0                   %% OPF_ALG not set, choose best option
    if have_fcn('pdipmopf')
      alg = 540;                %% PDIPM
    else
      alg = 560;                %% MIPS
    end
  end
  %% update deprecated algorithm codes to new, generalized equivalents
  if alg == 100 || alg == 200        %% CONSTR
    alg = 300;
  elseif alg == 120 || alg == 220    %% dense LP
    alg = 320;
  elseif alg == 140 || alg == 240    %% sparse (relaxed) LP
    alg = 340;
  elseif alg == 160 || alg == 260    %% sparse (full) LP
    alg = 360;
  end
  mpopt(11) = alg;
end

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

%% data dimensions
nb   = size(mpc.bus, 1);    %% number of buses
nl   = size(mpc.branch, 1); %% number of branches
ng   = size(mpc.gen, 1);    %% number of dispatchable injections
if isfield(mpc, 'A')
  nusr = size(mpc.A, 1);    %% number of linear user constraints
else
  nusr = 0;
end
if isfield(mpc, 'N')
  nw = size(mpc.N, 1);      %% number of general cost vars, w
else
  nw = 0;
end

%% add zero columns to bus, gen, branch for multipliers, etc if needed
if size(mpc.bus,2) < MU_VMIN
  mpc.bus = [mpc.bus zeros(nb, MU_VMIN-size(mpc.bus,2)) ];
end
if size(mpc.gen,2) < MU_QMIN
  mpc.gen = [ mpc.gen zeros(ng, MU_QMIN-size(mpc.gen,2)) ];
end
if size(mpc.branch,2) < MU_ANGMAX
  mpc.branch = [ mpc.branch zeros(nl, MU_ANGMAX-size(mpc.branch,2)) ];
end

if dc
  %% ignore reactive costs for DC
  mpc.gencost = pqcost(mpc.gencost, ng);

  %% reduce A and/or N from AC dimensions to DC dimensions, if needed
  if nusr || nw
    acc = [nb+(1:nb) 2*nb+ng+(1:ng)];   %% Vm and Qg columns
    if nusr && size(mpc.A, 2) >= 2*nb + 2*ng
      %% make sure there aren't any constraints on Vm or Qg
      if any(any(mpc.A(:, acc)))
        error('opf: attempting to solve DC OPF with user constraints on Vm or Qg');
      end
      mpc.A(:, acc) = [];               %% delete Vm and Qg columns
    end
    if nw && size(mpc.N, 2) >= 2*nb + 2*ng
      %% make sure there aren't any costs on Vm or Qg
      if any(any(mpc.N(:, acc)))
        error('opf: attempting to solve DC OPF with user costs on Vm or Qg');
      end
      mpc.N(:, acc) = [];               %% delete Vm and Qg columns
    end
  end
end

%% convert single-block piecewise-linear costs into linear polynomial cost
p1 = find(mpc.gencost(:, MODEL) == PW_LINEAR & mpc.gencost(:, NCOST) == 2);
% p1 = [];
if ~isempty(p1)
  x0 = mpc.gencost(p1, COST);
  y0 = mpc.gencost(p1, COST+1);
  x1 = mpc.gencost(p1, COST+2);
  y1 = mpc.gencost(p1, COST+3);
  m = (y1 - y0) ./ (x1 - x0);
  b = y0 - m .* x0;
  mpc.gencost(p1, MODEL) = POLYNOMIAL;
  mpc.gencost(p1, NCOST) = 2;
  mpc.gencost(p1, COST:COST+1) = [m b];
end

%% convert to internal numbering, remove out-of-service stuff
mpc = ext2int(mpc);

%% update dimensions
nb   = size(mpc.bus, 1);    %% number of buses
nl   = size(mpc.branch, 1); %% number of branches
ng   = size(mpc.gen, 1);    %% number of dispatchable injections

%% create (read-only) copies of individual fields for convenience
[baseMVA, bus, gen, branch, gencost, Au, lbu, ubu, mpopt, ...
    N, fparm, H, Cw, z0, zl, zu, userfcn] = opf_args(mpc, mpopt);

%% warn if there is more than one reference bus
refs = find(bus(:, BUS_TYPE) == REF);
if length(refs) > 1 && verbose > 0
  errstr = ['\nopf: Warning: Multiple reference buses.\n', ...
              '     For a system with islands, a reference bus in each island\n', ...
              '     may help convergence, but in a fully connected system such\n', ...
              '     a situation is probably not reasonable.\n\n' ];
  fprintf(errstr);
end

%% set up initial variables and bounds
Va   = bus(:, VA) * (pi/180);
Vm   = bus(:, VM);
Vm(gen(:, GEN_BUS)) = gen(:, VG);   %% buses with gens, init Vm from gen data
Pg   = gen(:, PG) / baseMVA;
Qg   = gen(:, QG) / baseMVA;
Pmin = gen(:, PMIN) / baseMVA;
Pmax = gen(:, PMAX) / baseMVA;
Qmin = gen(:, QMIN) / baseMVA;
Qmax = gen(:, QMAX) / baseMVA;

if dc               %% DC model
  %% more problem dimensions
  nv    = 0;            %% number of voltage magnitude vars
  nq    = 0;            %% number of Qg vars
  q1    = [];           %% index of 1st Qg column in Ay

  %% power mismatch constraints
  [B, Bf, Pbusinj, Pfinj] = makeBdc(baseMVA, bus, branch);
  neg_Cg = sparse(gen(:, GEN_BUS), 1:ng, -1, nb, ng);   %% Pbus w.r.t. Pg
  Amis = [B neg_Cg];
  bmis = -(bus(:, PD) + bus(:, GS)) / baseMVA - Pbusinj;

  %% branch flow constraints
  il = find(branch(:, RATE_A) ~= 0 & branch(:, RATE_A) < 1e10);
  nl2 = length(il);         %% number of constrained lines
  lpf = -Inf * ones(nl2, 1);
  upf = branch(il, RATE_A) / baseMVA - Pfinj(il);
  upt = branch(il, RATE_A) / baseMVA + Pfinj(il);

  user_vars = {'Va', 'Pg'};
  ycon_vars = {'Pg', 'y'};
else                %% AC model
  %% more problem dimensions
  nv    = nb;           %% number of voltage magnitude vars
  nq    = ng;           %% number of Qg vars
  q1    = 1+ng;         %% index of 1st Qg column in Ay

  %% dispatchable load, constant power factor constraints
  [Avl, lvl, uvl]  = makeAvl(baseMVA, gen);
  
  %% generator PQ capability curve constraints
  [Apqh, ubpqh, Apql, ubpql, Apqdata] = makeApq(baseMVA, gen);

  user_vars = {'Va', 'Vm', 'Pg', 'Qg'};
  ycon_vars = {'Pg', 'Qg', 'y'};
end

%% voltage angle reference constraints
Vau = Inf * ones(nb, 1);
Val = -Vau;
Vau(refs) = Va(refs);
Val(refs) = Va(refs);

%% branch voltage angle difference limits
[Aang, lang, uang, iang]  = makeAang(baseMVA, branch, nb, mpopt);

%% basin constraints for piece-wise linear gen cost variables
if alg == 545 || alg == 550     %% SC-PDIPM or TRALM, no CCV cost vars
  ny = 0;
  Ay = sparse(0, ng+nq);
  by =[];
else
  ipwl = find(gencost(:, MODEL) == PW_LINEAR);  %% piece-wise linear costs
  ny = size(ipwl, 1);   %% number of piece-wise linear cost vars
  [Ay, by] = makeAy(baseMVA, ng, gencost, 1, q1, 1+ng+nq);
end

%% more problem dimensions
nx    = nb+nv + ng+nq;  %% number of standard OPF control variables
if nusr
  nz = size(mpc.A, 2) - nx; %% number of user z variables
  if nz < 0
    error('opf: user supplied A matrix must have at least %d columns.', nx);
  end
else
  nz = 0;               %% number of user z variables
  if nw                 %% still need to check number of columns of N
    if size(mpc.N, 2) ~= nx;
      error('opf: user supplied N matrix must have %d columns.', nx);
    end
  end
end

%% construct OPF model object
om = opf_model(mpc);
if dc
  om = userdata(om, 'Bf', Bf);
  om = userdata(om, 'Pfinj', Pfinj);
  om = add_vars(om, 'Va', nb, Va, Val, Vau);
  om = add_vars(om, 'Pg', ng, Pg, Pmin, Pmax);
  om = add_constraints(om, 'Pmis', Amis, bmis, bmis, {'Va', 'Pg'}); %% nb
  om = add_constraints(om, 'Pf',  Bf(il,:), lpf, upf, {'Va'});      %% nl
  om = add_constraints(om, 'Pt', -Bf(il,:), lpf, upt, {'Va'});      %% nl
  om = add_constraints(om, 'ang', Aang, lang, uang, {'Va'});        %% nang
else
  om = add_vars(om, 'Va', nb, Va, Val, Vau);
  om = add_vars(om, 'Vm', nb, Vm, bus(:, VMIN), bus(:, VMAX));
  om = add_vars(om, 'Pg', ng, Pg, Pmin, Pmax);
  om = add_vars(om, 'Qg', ng, Qg, Qmin, Qmax);
  om = add_constraints(om, 'Pmis', nb, 'non-linear');
  om = add_constraints(om, 'Qmis', nb, 'non-linear');
  om = add_constraints(om, 'Sf', nl, 'non-linear');
  om = add_constraints(om, 'St', nl, 'non-linear');
  om = add_constraints(om, 'PQh', Apqh, [], ubpqh, {'Pg', 'Qg'});   %% npqh
  om = add_constraints(om, 'PQl', Apql, [], ubpql, {'Pg', 'Qg'});   %% npql
  om = add_constraints(om, 'vl',  Avl, lvl, uvl,   {'Pg', 'Qg'});   %% nvl
  om = add_constraints(om, 'ang', Aang, lang, uang, {'Va'});        %% nang
end

%% y vars, constraints for piece-wise linear gen costs
if ny > 0
  om = add_vars(om, 'y', ny);
  om = add_constraints(om, 'ycon', Ay, [], by, ycon_vars);          %% ncony
end

%% add user vars, constraints and costs (as specified via A, ..., N, ...)
if nz > 0
  om = add_vars(om, 'z', nz, z0, zl, zu);
  user_vars{end+1} = 'z';
end
if nusr
  om = add_constraints(om, 'usr', mpc.A, lbu, ubu, user_vars);      %% nusr
end
if nw
  user_cost.N = mpc.N;
  user_cost.Cw = Cw;
  if ~isempty(fparm)
    user_cost.dd = fparm(:, 1);
    user_cost.rh = fparm(:, 2);
    user_cost.kk = fparm(:, 3);
    user_cost.mm = fparm(:, 4);
  end
  if ~isempty(H)
    user_cost.H = H;
  end
  om = add_costs(om, 'usr', user_cost, user_vars);
end

%% execute userfcn callbacks for 'formulation' stage
om = run_userfcn(userfcn, 'formulation', om);

%% build user-defined costs
om = build_cost_params(om);

%% get indexing
[vv, ll, nn] = get_idx(om);

%% select optional solver output args
if nargout > 7 || nargout < 3
  output = struct('g', [], 'dg', []);
else
 output = struct([]);
end

%% call the specific solver
if verbose > 0
    v = mpver('all');
    fprintf('\nMATPOWER Version %s, %s', v.Version, v.Date);
end
if dc
  if verbose > 0
    fprintf(' -- DC Optimal Power Flow\n');
  end
  [results, success, raw] = dcopf_solver(om, mpopt, output);
else
  %%-----  call specific AC OPF solver  -----
  if verbose > 0
    fprintf(' -- AC Optimal Power Flow\n');
  end
  if alg == 500                                 %% MINOPF
    if ~have_fcn('minopf')
      error('opf: OPF_ALG %d requires MINOPF (see http://www.pserc.cornell.edu/minopf/)', alg);
    end
    [results, success, raw] = mopf_solver(om, mpopt, output);
  elseif alg == 300                             %% CONSTR
    if ~have_fcn('constr')
      error('opf: OPF_ALG %d requires CONSTR (Optimization Toolbox 1.x)', alg);
    end
    [results, success, raw] = copf_solver(om, mpopt, output);
  elseif alg == 320 || alg == 340 || alg == 360   %% LP
    [results, success, raw] = lpopf_solver(om, mpopt, output);
  elseif alg == 520                             %% FMINCON
    if ~have_fcn('fmincon')
      error('opf: OPF_ALG %d requires FMINCON (Optimization Toolbox 2.x or later)', alg);
    end
    if have_fcn('anon_fcns')
      solver = @fmincopf_solver;
    else
      solver = @fmincopf6_solver;
    end
    [results, success, raw] = feval(solver, om, mpopt, output);
  elseif alg == 540 || alg == 545 || alg == 550 %% PDIPM_OPF, SCPDIPM_OPF, or TRALM_OPF
    if alg == 540       % PDIPM_OPF
      if ~have_fcn('pdipmopf')
        error('opf: OPF_ALG %d requires PDIPMOPF (see http://www.pserc.cornell.edu/tspopf/)', alg);
      end
    elseif alg == 545       % SCPDIPM_OPF
      if ~have_fcn('scpdipmopf')
        error('opf: OPF_ALG %d requires SCPDIPMOPF (see http://www.pserc.cornell.edu/tspopf/)', alg);
      end
    elseif alg == 550       % TRALM_OPF
      if ~have_fcn('tralmopf')
        error('opf: OPF_ALG %d requires TRALM (see http://www.pserc.cornell.edu/tspopf/)', alg);
      end
    end
    [results, success, raw] = tspopf_solver(om, mpopt, output);
  elseif alg == 560 || alg == 565               %% MIPS
    if have_fcn('anon_fcns')
      solver = @mipsopf_solver;
    else
      solver = @mips6opf_solver;
    end
    [results, success, raw] = feval(solver, om, mpopt, output);
  end
end
if ~isfield(raw, 'output') || ~isfield(raw.output, 'alg') || isempty(raw.output.alg)
    raw.output.alg = alg;
end
if isfield(results, 'g')
  g = results.g;
end
if isfield(results, 'dg')
  jac = results.dg;
end

if success
  if ~dc
    %% copy bus voltages back to gen matrix
    results.gen(:, VG) = results.bus(results.gen(:, GEN_BUS), VM);
  
    %% gen PQ capability curve multipliers
    if ll.N.PQh > 0 || ll.N.PQl > 0
      mu_PQh = results.mu.lin.l(ll.i1.PQh:ll.iN.PQh) - results.mu.lin.u(ll.i1.PQh:ll.iN.PQh);
      mu_PQl = results.mu.lin.l(ll.i1.PQl:ll.iN.PQl) - results.mu.lin.u(ll.i1.PQl:ll.iN.PQl);
      results.gen = update_mupq(baseMVA, results.gen, mu_PQh, mu_PQl, Apqdata);
    end
  end

  %% angle limit constraint multipliers
  if ll.N.ang > 0
    results.branch(iang, MU_ANGMIN) = results.mu.lin.l(ll.i1.ang:ll.iN.ang) * pi/180;
    results.branch(iang, MU_ANGMAX) = results.mu.lin.u(ll.i1.ang:ll.iN.ang) * pi/180;
  end
end

%% assign values and limit shadow prices for variables
om_var_order = get(om, 'var', 'order');
for k = 1:length(om_var_order)
  name = om_var_order{k};
  if getN(om, 'var', name)
    idx = vv.i1.(name):vv.iN.(name);
    results.var.val.(name) = results.x(idx);
    results.var.mu.l.(name) = results.mu.var.l(idx);
    results.var.mu.u.(name) = results.mu.var.u(idx);
  end
end

%% assign shadow prices for linear constraints
om_lin_order = get(om, 'lin', 'order');
for k = 1:length(om_lin_order)
  name = om_lin_order{k};
  if getN(om, 'lin', name)
    idx = ll.i1.(name):ll.iN.(name);
    results.lin.mu.l.(name) = results.mu.lin.l(idx);
    results.lin.mu.u.(name) = results.mu.lin.u(idx);
  end
end

%% assign shadow prices for non-linear constraints
if ~dc
  om_nln_order = get(om, 'nln', 'order');
  for k = 1:length(om_nln_order)
    name = om_nln_order{k};
    if getN(om, 'nln', name)
      idx = nn.i1.(name):nn.iN.(name);
      results.nln.mu.l.(name) = results.mu.nln.l(idx);
      results.nln.mu.u.(name) = results.mu.nln.u(idx);
    end
  end
end

%% assign values for components of user cost
om_cost_order = get(om, 'cost', 'order');
for k = 1:length(om_cost_order)
  name = om_cost_order{k};
  if getN(om, 'cost', name)
    results.cost.(name) = compute_cost(om, results.x, name);
  end
end

%% revert to original ordering, including out-of-service stuff
results = int2ext(results);

%% zero out result fields of out-of-service gens & branches
if ~isempty(results.order.gen.status.off)
  results.gen(results.order.gen.status.off, [PG QG MU_PMAX MU_PMIN]) = 0;
end
if ~isempty(results.order.branch.status.off)
  results.branch(results.order.branch.status.off, [PF QF PT QT MU_SF MU_ST MU_ANGMIN MU_ANGMAX]) = 0;
end

%% if single-block PWL costs were converted to POLY, insert dummy y into x
%% Note: The "y" portion of x will be messed up, but everything else
%%       should be in the expected locations.
if ~isempty(p1) && alg ~= 545 && alg ~= 550
  if dc
    nx = vv.N.Pg;
  else
    nx = vv.N.Qg;
  end
  y = zeros(length(p1), 1);
  raw.xr = [ raw.xr(1:nx); y; raw.xr(nx+1:end)];
  results.x = [ results.x(1:nx); y; results.x(nx+1:end)];
end

%% compute elapsed time
et = etime(clock, t0);

%% finish preparing output
if nargout > 0
  if nargout <= 2
    results.et = et;
    results.success = success;
    results.raw = raw;
    busout = results;
    genout = success;
  else
    [busout, genout, branchout, f, info, xr, pimul] = deal(results.bus, ...
        results.gen, results.branch, results.f, raw.info, raw.xr, raw.pimul);
  end
elseif success
  results.et = et;
  results.success = success;
  printpf(results, 1, mpopt);
end

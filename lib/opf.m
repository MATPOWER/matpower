function [busout, genout, branchout, f, success, info, et, g, jac, xr, pimul] = ...
    opf(varargin)
%OPF  Solves an optimal power flow.
%
%   For an AC OPF, if the OPF algorithm is not set explicitly in the options,
%   it will choose the best available solver, searching in the following order:
%   minopf, pdipmopf, fmincon, LP-base, and constr.
%
%   Returns either a results struct and an optional success flag, or individual
%   data matrices, the objective function value and a success flag. In the
%   latter case, there are additional optional return values.
%
%   results = opf(...)
%   [results, success] = opf(...)
%   [bus, gen, branch, f, success] = opf(...)
%   [bus, gen, branch, f, success, info, et, g, jac, xr, pimul] = opf(...)
%
%   Input arguments options are as follows:
%
%   opf(mpc)
%   opf(mpc, mpopt)
%   opf(mpc, userfcn, mpopt)
%   opf(mpc, A, l, u)
%   opf(mpc, A, l, u, mpopt)
%   opf(mpc, A, l, u, mpopt, N, fparm, H, Cw)
%   opf(mpc, A, l, u, mpopt, N, fparm, H, Cw, z0, zl, zu)
%
%   opf(baseMVA, bus, gen, branch, areas, gencost)
%   opf(baseMVA, bus, gen, branch, areas, gencost, mpopt)
%   opf(baseMVA, bus, gen, branch, areas, gencost, userfcn, mpopt)
%   opf(baseMVA, bus, gen, branch, areas, gencost, A, l, u)
%   opf(baseMVA, bus, gen, branch, areas, gencost, A, l, u, mpopt)
%   opf(baseMVA, bus, gen, branch, areas, gencost, A, l, u, ...
%                               mpopt, N, fparm, H, Cw)
%   opf(baseMVA, bus, gen, branch, areas, gencost, A, l, u, ...
%                               mpopt, N, fparm, H, Cw, z0, zl, zu)
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

%   MATPOWER
%   $Id$
%   by Ray Zimmerman, PSERC Cornell
%   and Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Autonoma de Manizales
%   Copyright (c) 1996-2008 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%%----- initialization -----
t0 = clock;         %% start timer

%% process input arguments
[baseMVA, bus, gen, branch, areas, gencost, Au, lbu, ubu, mpopt, ...
    N, fparm, H, Cw, z0, zl, zu, userfcn] = opf_args(varargin{:});

%% options
dc  = mpopt(10);        %% PF_DC        : 1 = DC OPF, 0 = AC OPF 
alg = mpopt(11);        %% OPF_ALG

%% set AC OPF algorithm code
if ~dc
  if alg == 0                   %% OPF_ALG not set, choose best option
    if have_fcn('minopf')
      alg = 500;                %% MINOS
    elseif have_fcn('pdipmopf')
      alg = 540;                %% PDIPM
    else
      alg = 560;
    else
      error('opf: no OPF solvers available');
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
nb   = size(bus, 1);    %% number of buses
nl   = size(branch, 1); %% number of branches
ng   = size(gen, 1);    %% number of dispatchable injections
nusr = size(Au, 1);     %% number of linear user constraints
nw   = size(N, 1);      %% number of general cost vars, w

%% add zero columns to bus, gen, branch for multipliers, etc if needed
if size(bus,2) < MU_VMIN
  bus = [bus zeros(nb ,MU_VMIN-size(bus,2)) ];
end
if size(gen,2) < MU_QMIN
  gen = [ gen zeros(ng,MU_QMIN-size(gen,2)) ];
end
if size(branch,2) < MU_ANGMAX
  branch = [ branch zeros(nl,MU_ANGMAX-size(branch,2)) ];
end

if dc
  %% ignore reactive costs for DC
  gencost = pqcost(gencost, ng);

  %% reduce Au and/or N from AC dimensions to DC dimensions, if needed
  if nusr || nw
    acc = [nb+(1:nb) 2*nb+ng+(1:ng)];   %% Vm and Qg columns
    if nusr && size(Au, 2) >= 2*nb + 2*ng
      %% make sure there aren't any constraints on Vm or Qg
      if any(any(Au(:, acc)))
        error('opf: attempting to solve DC OPF with user constraints on Vm or Qg');
      end
      Au(:, acc) = [];                  %% delete Vm and Qg columns
    end
    if nw && size(N, 2) >= 2*nb + 2*ng
      %% make sure there aren't any costs on Vm or Qg
      if any(any(N(:, acc)))
        error('opf: attempting to solve DC OPF with user constraints on Vm or Qg');
      end
      N(:, acc) = [];                   %% delete Vm and Qg columns
    end
  end
end

%% filter out inactive generators and branches; save original bus & branch
status.gen.on     = find(gen(:,GEN_STATUS) > 0);
status.gen.off    = find(gen(:,GEN_STATUS) <= 0);
status.branch.on  = find(branch(:,BR_STATUS) ~= 0);
status.branch.off = find(branch(:,BR_STATUS) == 0);
original.gen      = gen;            %% save originals
original.branch   = branch;
gen(      status.gen.off, :) = [];      %% delete out-of-service gens
branch(status.branch.off, :) = [];      %% delete out-of-service branches
if size(gencost,1) == ng                %% delete costs for out-of-service gens
  gencost(status.gen.off, :) = [];
else
  gencost([status.gen.off; status.gen.off+ng], :) = [];
end
if dc
  offcols = status.gen.off + nb;
else
  offcols = [status.gen.off; status.gen.off+ng] + 2*nb;
end
if nusr
  if any(any(Au(:, offcols)))
    error('opf: user constraint involves out-of-service gen');
  end
  Au(:, offcols) = [];
end
if nw
  N(:, offcols) = [];
end

%% update dimensions
ng = size(gen, 1);      %% number of dispatchable injections
nl = size(branch, 1);   %% number of branches

%% convert to internal consecutive bus numbering
[reorder.i2e, bus, gen, branch] = ext2int(bus, gen, branch);

%% sort generators in order of increasing bus number
[tmp, igen] = sort(gen(:, GEN_BUS));
[tmp, reorder.invgen] = sort(igen);    %% save for inverse reordering later
reorder.gen = igen;
gen  = gen(igen, :);
if ng == size(gencost,1)
  gencost = gencost(igen, :);
else
  gencost = gencost( [igen; igen+ng], :);
end
one2ng = (1:ng)';
if dc
  oldcols = one2ng + nb;
  newcols =   igen + nb;
else
  oldcols = [one2ng; one2ng+ng] + 2*nb;
  newcols = [  igen;   igen+ng] + 2*nb;
end
if nusr
  Au(:, oldcols) = Au(:, newcols);
end
if nw
  N(:, oldcols) = N(:, newcols);
end

%% warn if there is more than one reference bus
refs = find(bus(:, BUS_TYPE) == REF);
if length(refs) > 1
  errstr = ['\nopf: Warning: Multiple reference buses.\n', ...
              '     For a system with islands, a reference bus in each island\n', ...
              '     may help convergence, but in a fully connected system such\n', ...
              '     a situation is probably not reasonable.\n\n' ];
  fprintf(errstr);
end

%% set up MATPOWER case data struct
mpc = struct('baseMVA', baseMVA, 'bus', bus, 'gen', gen, ...
    'branch', branch, 'gencost', gencost, 'original', original, ...
    'status', status, 'reorder', reorder);
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
  mpc.Bf = Bf;
  mpc.Pfinj = Pfinj;
  neg_Cg = sparse(gen(:, GEN_BUS), 1:ng, -1, nb, ng);   %% Pbus w.r.t. Pg
  Amis = [B neg_Cg];
  bmis = -(bus(:, PD) + bus(:, GS)) / baseMVA - Pbusinj;

  %% branch flow constraints
  il = find(branch(:, RATE_A) ~= 0 & branch(:, RATE_A) < 1e10);
  nl2 = length(il);			%% number of constrained lines
  lpf = -Inf * ones(nl2, 1);
  upf = branch(il, RATE_A) / baseMVA - Pfinj(il);
  upt = branch(il, RATE_A) / baseMVA + Pfinj(il);

  %% voltage angel reference constraints
  Vau = Inf * ones(nb, 1);
  Val = -Vau;
  Vau(refs) = Va(refs);
  Val(refs) = Va(refs);
  
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
if isempty(Au)      %% set nz
  nz = 0;               %% number of user z variables
  Au = sparse(0,nx);
  if ~isempty(N)        %% still need to check number of columns of N
    if size(N, 2) ~= nx;
      error('opf: user supplied N matrix must have %d columns.', nx);
    end
  end
else
  nz = size(Au,2) - nx; %% number of user z variables
  if nz < 0
    error('opf: user supplied A matrix must have at least %d columns.', nx);
  end
end

%% set up user cost params
user_cost.N = N;
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

%% construct OPF model object
om = opf_model(mpc);
if dc
  om = add_vars(om, 'Va', nb, Va, Val, Vau);
  om = add_vars(om, 'Pg', ng, Pg, Pmin, Pmax);
  om = add_constraints(om, 'Pmis', Amis, bmis, bmis, {'Va', 'Pg'}); %% nb
  om = add_constraints(om, 'Pf',  Bf(il,:), lpf, upf, {'Va'});      %% nl
  om = add_constraints(om, 'Pt', -Bf(il,:), lpf, upt, {'Va'});      %% nl
  om = add_constraints(om, 'ang', Aang, lang, uang, {'Va'});        %% nang
else
  om = add_vars(om, 'Va', nb, Va);
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

%% add user vars, constraints and costs (as specified via Au, ..., N, ...)
if nz > 0
  om = add_vars(om, 'z', nz, z0, zl, zu);
  user_vars{end+1} = 'z';
end
if nusr
  om = add_constraints(om, 'usr', Au, lbu, ubu, user_vars);         %% nusr
end
if nw
  om = add_costs(om, 'usr', user_cost, user_vars);
end

%% add user vars, constraints, costs (as specified via userfcn)
if ~isempty(userfcn) && isfield(userfcn, 'name')
  for k = 1:length(userfcn)
	if isfield(userfcn, 'args') && ~isempty(userfcn(k).args)
	  om = feval(userfcn(k).name, om, userfcn(k).args);
	else
	  om = feval(userfcn(k).name, om);
	end
  end
end

%% get indexing
[vv, ll, nn] = get_idx(om);

%% select optional solver output args
if nargout > 7
  output = struct('g', [], 'dg', []);
else
 output = struct([]);
end

%% call the specific solver
if dc
  [results, success, raw] = dcopf_solver(om, mpopt, output);
%   pimul = [ ...
%       results.mu.lin.l - results.mu.lin.u;
%       -ones(ny>0, 1);
%       results.mu.var.l - results.mu.var.u;
%   ];
else
  %%-----  call specific AC OPF solver  -----
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
    mlver = ver('matlab');
    if str2double(mlver.Version(1)) < 7    %% anonymous functions not available
      fmc = @fmincopf6_solver;
    else
      fmc = @fmincopf_solver;
    end
    [results, success, raw] = feval(fmc, om, mpopt, output);
  elseif alg == 540 || alg == 545 || alg == 550   %% PDIPM_OPF, SCPDIPM_OPF, or TRALM_OPF
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
  elseif alg == 560                             %% pdipm (pure Matlab)
    [results, success, raw] = pdipm_solver(om, mpopt, output);
  end
%   pimul = [ ...
%       results.mu.nln.l - results.mu.nln.u;
%       results.mu.lin.l - results.mu.lin.u;
%       -ones(ny>0, 1);
%       results.mu.var.l - results.mu.var.u;
%   ];
end
[bus, gen, branch, f, info, xr, pimul] = deal(results.bus, results.gen, ...
                results.branch, results.f, raw.info, raw.xr, raw.pimul);
if isfield(results, 'g')
  g = results.g;
end
if isfield(results, 'dg')
  jac = results.dg;
end
% xr = results.x;
xr = raw.xr;

% norm(xr - raw.xr(1:length(xr)))
% norm(pimul - raw.pimul(1:length(pimul)))
% fprintf('%g\t%g\t%g\n', [pimul raw.pimul abs(pimul - raw.pimul)]');

%% gen PQ capability curve multipliers
if ~dc && success && (ll.N.PQh > 0 || ll.N.PQl > 0)
  mu_PQh = results.mu.lin.l(ll.i1.PQh:ll.iN.PQh) - results.mu.lin.u(ll.i1.PQh:ll.iN.PQh);
  mu_PQl = results.mu.lin.l(ll.i1.PQl:ll.iN.PQl) - results.mu.lin.u(ll.i1.PQl:ll.iN.PQl);
  gen = update_mupq(baseMVA, gen, mu_PQh, mu_PQl, Apqdata);
end

%% angle limit constraint multipliers
if success && (ll.N.ang > 0)
  branch(iang, MU_ANGMIN) = results.mu.lin.l(ll.i1.ang:ll.iN.ang) * pi/180;
  branch(iang, MU_ANGMAX) = results.mu.lin.u(ll.i1.ang:ll.iN.ang) * pi/180;
end

%% revert to original gen ordering
gen = gen(reorder.invgen, :);
gen(:, VG) = bus(gen(:, GEN_BUS), VM);  %% copy bus voltages back to gen matrix

%% convert to original external bus ordering
[bus, gen, branch] = int2ext(reorder.i2e, bus, gen, branch);

%% include out-of-service branches and gens
genout    = original.gen;
branchout = original.branch;
genout(status.gen.on, :) = gen;
branchout(status.branch.on, :)  = branch;
if ~isempty(status.gen.off)     %% zero out result fields of out-of-service gens
  genout(status.gen.off, [PG QG MU_PMAX MU_PMIN]) = 0;
end
if ~isempty(status.branch.off)  %% zero out result fields of out-of-service branches
  branchout(status.branch.off, [PF QF PT QT MU_SF MU_ST MU_ANGMIN MU_ANGMAX]) = 0;
end

%% compute elapsed time
et = etime(clock, t0);

%% finish preparing output
if nargout > 0
  if nargout <= 2
    results.bus = bus;
    results.gen = genout;
    results.branch = branchout;
    results.status = status;
    results.reorder = reorder;
    results.et = et;
    results.raw = raw;
    results.success = success;
    results.om = om;
    
    %% values and limit shadow prices for variables
    om_var_order = get(om, 'var', 'order');
    for k = 1:length(om_var_order)
      name = om_var_order{k};
      if get_var_N(om, name)
        idx = vv.i1.(name):vv.iN.(name);
        results.var.val.(name) = results.x(idx);
        results.var.mu.l.(name) = results.mu.var.l(idx);
        results.var.mu.u.(name) = results.mu.var.u(idx);
      end
    end

    %% shadow prices for linear constraints
    om_lin_order = get(om, 'lin', 'order');
    for k = 1:length(om_lin_order)
      name = om_lin_order{k};
      if get_lin_N(om, name)
        idx = ll.i1.(name):ll.iN.(name);
        results.lin.mu.l.(name) = results.mu.lin.l(idx);
        results.lin.mu.u.(name) = results.mu.lin.u(idx);
      end
    end

    %% shadow prices for non-linear constraints
    if ~dc
      om_nln_order = get(om, 'nln', 'order');
      for k = 1:length(om_nln_order)
        name = om_nln_order{k};
        if get_nln_N(om, name)
          idx = nn.i1.(name):nn.iN.(name);
          results.nln.mu.l.(name) = results.mu.nln.l(idx);
          results.nln.mu.u.(name) = results.mu.nln.u(idx);
        end
      end
    end

    busout = results;
    genout = success;
  else
    busout = bus;
  end
elseif success
  printpf(baseMVA, bus, genout, branchout, f, success, et, 1, mpopt);
end

return;

function om = opf_setup(mpc, mpopt)
%OPF  Constructs an OPF model object from a MATPOWER case struct.
%   OM = OPF_SETUP(MPC, MPOPT)
%
%   Assumes that MPC is a MATPOWER case struct with internal indexing,
%   all equipment in-service, etc.
%
%   See also OPF, EXT2INT, OPF_EXECUTE.

%   MATPOWER
%   Copyright (c) 1996-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%   and Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Nacional de Colombia
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% options
dc  = strcmp(upper(mpopt.model), 'DC');
alg = upper(mpopt.opf.ac.solver);
use_vg = mpopt.opf.use_vg;
vcart = ~dc && mpopt.opf.v_cartesian;

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

%% define flag to indicate whether we are tied to legacy formulation
%% implemented by legacy MINOS and PDIPM solvers (e.g. with hard-coded
%% costs and constrain
if strcmp(alg, 'MINOPF') || strcmp(alg, 'PDIPM') || ...
        strcmp(alg, 'TRALM') || strcmp(alg, 'SDPOPF')
    legacy_formulation = 1;
    if vcart
        error('Option ''opf.v_cartesian'' = 1 is not compatible with ''opf.solver.ac''=''%s''.', alg);
    end
    if mpopt.opf.current_balance
        error('Option ''opf.current_balance'' = 1 is not compatible with ''opf.solver.ac''=''%s''.', alg);
    end
else
    legacy_formulation = 0;
end
if ~dc && ( ~isempty(mpopt.exp.sys_wide_zip_loads.pw) && ...
                    ~isequal(mpopt.exp.sys_wide_zip_loads.pw, [1 0 0]) || ...
            ~isempty(mpopt.exp.sys_wide_zip_loads.qw) && ...
                    ~isequal(mpopt.exp.sys_wide_zip_loads.qw, [1 0 0]) )
    if vcart
        warning('Voltage dependent loads are not supported with option ''opf.v_cartesian'' = 1. Reverting to constant power load model.');
    end
    if mpopt.opf.current_balance
        warning('Voltage dependent loads are not supported with option ''opf.current_balance'' = 1. Reverting to constant power load model.');
    end
end

%% data dimensions
nb   = size(mpc.bus, 1);    %% number of buses
nl   = size(mpc.branch, 1); %% number of branches
ng   = size(mpc.gen, 1);    %% number of dispatchable injections
nnle = 0;                   %% number of nonlinear user-defined equality cons
nnli = 0;                   %% number of nonlinear user-defined inequality cons
if isfield(mpc, 'A')
  nlin = size(mpc.A, 1);    %% number of linear user constraints
else
  nlin = 0;
end
if isfield(mpc, 'N')
  nw = size(mpc.N, 1);      %% number of general cost vars, w
else
  nw = 0;
end

if dc
  %% ignore reactive costs for DC
  mpc.gencost = pqcost(mpc.gencost, ng);

  %% reduce A and/or N from AC dimensions to DC dimensions, if needed
  if nlin || nw
    acc = [nb+(1:nb) 2*nb+ng+(1:ng)];   %% Vm and Qg columns
    if nlin && size(mpc.A, 2) >= 2*nb + 2*ng
      %% make sure there aren't any constraints on Vm or Qg
      if any(any(mpc.A(:, acc)))
        error('opf_setup: attempting to solve DC OPF with user constraints on Vm or Qg');
      end
      mpc.A(:, acc) = [];               %% delete Vm and Qg columns
    end
    if nw && size(mpc.N, 2) >= 2*nb + 2*ng
      %% make sure there aren't any costs on Vm or Qg
      if any(any(mpc.N(:, acc)))
        [ii, jj] = find(mpc.N(:, acc));
        ii = unique(ii);    %% indices of w with potential non-zero cost terms from Vm or Qg
        if any(mpc.Cw(ii)) || (isfield(mpc, 'H') && ~isempty(mpc.H) && ...
                any(any(mpc.H(:, ii))))
          error('opf_setup: attempting to solve DC OPF with user costs on Vm or Qg');
        end
      end
      mpc.N(:, acc) = [];               %% delete Vm and Qg columns
    end
  end
else    %% AC
  if use_vg     %% adjust bus voltage limits based on generator Vg setpoint
    %% gen connection matrix, element i, j is 1 if, generator j at bus i is ON
    Cg = sparse(mpc.gen(:, GEN_BUS), (1:ng)', mpc.gen(:, GEN_STATUS) > 0, nb, ng);
    Vbg = Cg * sparse(1:ng, 1:ng, mpc.gen(:, VG), ng, ng);
    Vmax = max(Vbg, [], 2); %% zero for non-gen buses, else max Vg of gens @ bus
    ib = find(Vmax);                %% buses with online gens
    Vmin = max(2*Cg - Vbg, [], 2);  %% same as Vmax, except min Vg of gens @ bus
    Vmin(ib) = 2 - Vmin(ib);

    if use_vg == 1      %% use Vg setpoint directly
        mpc.bus(ib, VMAX) = Vmax(ib);   %% max set by max Vg @ bus
        mpc.bus(ib, VMIN) = Vmin(ib);   %% min set by min Vg @ bus
        mpc.bus(ib, VM) = mpc.bus(ib, VMAX);
    elseif use_vg > 0 && use_vg < 1     %% fractional value
        %% use weighted avg between original Vmin/Vmax limits and Vg
        mpc.bus(ib, VMAX) = (1-use_vg) * mpc.bus(ib, VMAX) + use_vg * Vmax(ib);
        mpc.bus(ib, VMIN) = (1-use_vg) * mpc.bus(ib, VMIN) + use_vg * Vmin(ib);
    else
        error('opf_setup: option ''opf.use_vg'' (= %g) cannot be negative or greater than 1', use_vg);
    end
  end
  if isfield(mpc, 'user_constraints')
    if isfield(mpc.user_constraints, 'nle')
      for k = 1:length(mpc.user_constraints.nle)
        nnle = nnle + mpc.user_constraints.nle{k}{2};
      end
    end
    if isfield(mpc.user_constraints, 'nli')
      for k = 1:length(mpc.user_constraints.nli)
        nnli = nnli + mpc.user_constraints.nli{k}{2};
      end
    end
  end
end

%% convert single-block piecewise-linear costs into linear polynomial cost
pwl1 = find(mpc.gencost(:, MODEL) == PW_LINEAR & mpc.gencost(:, NCOST) == 2);
% p1 = [];
if ~isempty(pwl1)
  x0 = mpc.gencost(pwl1, COST);
  y0 = mpc.gencost(pwl1, COST+1);
  x1 = mpc.gencost(pwl1, COST+2);
  y1 = mpc.gencost(pwl1, COST+3);
  m = (y1 - y0) ./ (x1 - x0);
  b = y0 - m .* x0;
  mpc.gencost(pwl1, MODEL) = POLYNOMIAL;
  mpc.gencost(pwl1, NCOST) = 2;
  mpc.gencost(pwl1, COST:COST+1) = [m b];
end

%% create (read-only) copies of individual fields for convenience
[baseMVA, bus, gen, branch, gencost, Au, lbu, ubu, mpopt, ...
    N, fparm, H, Cw, z0, zl, zu, userfcn] = opf_args(mpc, mpopt);

%% warn if there is more than one reference bus
refs = find(bus(:, BUS_TYPE) == REF);
if length(refs) > 1 && mpopt.verbose > 0
  errstr = ['\nopf_setup: Warning: Multiple reference buses.\n', ...
              '           For a system with islands, a reference bus in each island\n', ...
              '           may help convergence, but in a fully connected system such\n', ...
              '           a situation is probably not reasonable.\n\n' ];
  fprintf(errstr);
end

%% set up initial variables and bounds
Va   = bus(:, VA) * (pi/180);
Vau = Inf(nb, 1);       %% voltage angle limits
Val = -Vau;
Vau(refs) = Va(refs);   %% voltage angle reference constraints
Val(refs) = Va(refs);
Pg   = gen(:, PG) / baseMVA;
Pmin = gen(:, PMIN) / baseMVA;
Pmax = gen(:, PMAX) / baseMVA;
if ~dc
  Vm   = bus(:, VM);
  Qg   = gen(:, QG) / baseMVA;
  Qmin = gen(:, QMIN) / baseMVA;
  Qmax = gen(:, QMAX) / baseMVA;
  if vcart
    V = Vm .* exp(1j*Va);
    Vr = real(V);
    Vi = imag(V);
  end
end

%% find/prepare polynomial generator costs
cpg = [];
cqg = [];
[pcost qcost] = pqcost(mpc.gencost, ng);
ip0 = find(pcost(:, MODEL) == POLYNOMIAL & pcost(:, NCOST) == 1);   %% constant
ip1 = find(pcost(:, MODEL) == POLYNOMIAL & pcost(:, NCOST) == 2);   %% linear
ip2 = find(pcost(:, MODEL) == POLYNOMIAL & pcost(:, NCOST) == 3);   %% quadratic
ip3 = find(pcost(:, MODEL) == POLYNOMIAL & pcost(:, NCOST) > 3);    %% cubic or greater
if ~isempty(ip2) || ~isempty(ip1) || ~isempty(ip0)
    kpg = zeros(ng, 1);
    cpg = zeros(ng, 1);
    if ~isempty(ip2)
        Qpg = zeros(ng, 1);
        Qpg(ip2) = 2 * pcost(ip2, COST) * baseMVA^2;
        cpg(ip2) = cpg(ip2) + pcost(ip2, COST+1) * baseMVA;
        kpg(ip2) = kpg(ip2) + pcost(ip2, COST+2);
    else
        Qpg = [];   %% no quadratic terms
    end
    cpg(ip1) = cpg(ip1) + pcost(ip1, COST) * baseMVA;
    kpg(ip1) = kpg(ip1) + pcost(ip1, COST+1);
    kpg(ip0) = kpg(ip0) + pcost(ip0, COST);
end
if ~isempty(qcost)
    iq0 = find(qcost(:, MODEL) == POLYNOMIAL & qcost(:, NCOST) == 1);   %% constant
    iq1 = find(qcost(:, MODEL) == POLYNOMIAL & qcost(:, NCOST) == 2);   %% linear
    iq2 = find(qcost(:, MODEL) == POLYNOMIAL & qcost(:, NCOST) == 3);   %% quadratic
    iq3 = find(qcost(:, MODEL) == POLYNOMIAL & qcost(:, NCOST) > 3);    %% cubic or greater
    if ~isempty(iq2) || ~isempty(iq1) || ~isempty(iq0)
        kqg = zeros(ng, 1);
        cqg = zeros(ng, 1);
        if ~isempty(iq2)
            Qqg = zeros(ng, 1);
            Qqg(iq2) = 2 * qcost(iq2, COST) * baseMVA^2;
            cqg(iq2) = cqg(iq2) + qcost(iq2, COST+1) * baseMVA;
            kqg(iq2) = kqg(iq2) + qcost(iq2, COST+2);
        else
            Qqg = [];   %% no quadratic terms
        end
        cqg(iq1) = cqg(iq1) + qcost(iq1, COST) * baseMVA;
        kqg(iq1) = kqg(iq1) + qcost(iq1, COST+1);
        kqg(iq0) = kqg(iq0) + qcost(iq0, COST);
    end
end

%% branch voltage angle difference limits
[Aang, lang, uang, iang]  = makeAang(baseMVA, branch, nb, mpopt);
nang = length(iang);

if dc               %% DC model
  %% check generator costs
  if ~isempty(ip3)
    error('opf_setup: DC OPF cannot handle polynomial costs with higher than quadratic order.');
  end

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
  if nl2
    upf = branch(il, RATE_A) / baseMVA - Pfinj(il);
    upt = branch(il, RATE_A) / baseMVA + Pfinj(il);
  else
    upf = [];
    upt = [];
  end

  user_vars = {'Va', 'Pg'};
  ycon_vars = {'Pg', 'y'};
else                %% AC model
  %% more problem dimensions
  nv    = nb;           %% number of voltage magnitude vars
  nq    = ng;           %% number of Qg vars
  q1    = 1+ng;         %% index of 1st Qg column in Ay

  %% find branches with flow limits
  il = find(branch(:, RATE_A) ~= 0 & branch(:, RATE_A) < 1e10);
  nl2 = length(il);         %% number of constrained lines

  %% build admittance matrices
  [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);

  %% dispatchable load, constant power factor constraints
  [Avl, lvl, uvl]  = makeAvl(mpc);

  %% generator PQ capability curve constraints
  [Apqh, ubpqh, Apql, ubpql, Apqdata] = makeApq(baseMVA, gen);

  if vcart
      user_vars = {'Vr', 'Vi', 'Pg', 'Qg'};
      nodal_balance_vars = {'Vr', 'Vi', 'Pg', 'Qg'};
      flow_lim_vars = {'Vr', 'Vi'};
  else
      user_vars = {'Va', 'Vm', 'Pg', 'Qg'};
      nodal_balance_vars = {'Va', 'Vm', 'Pg', 'Qg'};
      flow_lim_vars = {'Va', 'Vm'};
  end
  ycon_vars = {'Pg', 'Qg', 'y'};

  %% nonlinear constraint functions
  if mpopt.opf.current_balance
    mis_cons = {'rImis', 'iImis'};
    fcn_mis = @(x)opf_current_balance_fcn(x, mpc, Ybus, mpopt);
    hess_mis = @(x, lam)opf_current_balance_hess(x, lam, mpc, Ybus, mpopt);
  else
    mis_cons = {'Pmis', 'Qmis'};
    fcn_mis = @(x)opf_power_balance_fcn(x, mpc, Ybus, mpopt);
    hess_mis = @(x, lam)opf_power_balance_hess(x, lam, mpc, Ybus, mpopt);
  end
  fcn_flow = @(x)opf_branch_flow_fcn(x, mpc, Yf(il, :), Yt(il, :), il, mpopt);
  hess_flow = @(x, lam)opf_branch_flow_hess(x, lam, mpc, Yf(il, :), Yt(il, :), il, mpopt);
  if vcart
    fcn_vref = @(x)opf_vref_fcn(x, mpc, refs, mpopt);
    hess_vref = @(x, lam)opf_vref_hess(x, lam, mpc, refs, mpopt);
    veq = find(mpc.bus(:, VMIN) == mpc.bus(:, VMAX));
    viq = find(mpc.bus(:, VMIN) ~= mpc.bus(:, VMAX));
    nveq = length(veq);
    nvlims = length(viq);
    if nveq
      fcn_veq = @(x)opf_veq_fcn(x, mpc, veq, mpopt);
      hess_veq = @(x, lam)opf_veq_hess(x, lam, mpc, veq, mpopt);
    end
    fcn_vlim = @(x)opf_vlim_fcn(x, mpc, viq, mpopt);
    hess_vlim = @(x, lam)opf_vlim_hess(x, lam, mpc, viq, mpopt);
    fcn_ang = @(x)opf_branch_ang_fcn(x, Aang, lang, uang, iang, mpopt);
    hess_ang = @(x, lam)opf_branch_ang_hess(x, lam, Aang, lang, uang, iang, mpopt);
  end
  
  %% nonlinear cost functions
  if ~isempty(ip3)
    cost_Pg = @(x)opf_gen_cost_fcn(x, baseMVA, pcost, ip3, mpopt);
  end
  if ~isempty(qcost) && ~isempty(iq3)
    cost_Qg = @(x)opf_gen_cost_fcn(x, baseMVA, qcost, iq3, mpopt);
  end
end

%% basin constraints for piece-wise linear gen cost variables
if (strcmp(alg, 'PDIPM') && mpopt.pdipm.step_control) || strcmp(alg, 'TRALM')
  %% SC-PDIPM or TRALM, no CCV cost vars
  ny = 0;
  Ay = sparse(0, ng+nq);
  by =[];
else
  ipwl = find(gencost(:, MODEL) == PW_LINEAR);  %% piece-wise linear costs
  ny = size(ipwl, 1);   %% number of piece-wise linear cost vars
  [Ay, by] = makeAy(baseMVA, ng, gencost, 1, q1, 1+ng+nq);
end
if any(gencost(:, MODEL) ~= POLYNOMIAL & gencost(:, MODEL) ~= PW_LINEAR)
    error('opf_setup: some generator cost rows have invalid MODEL value');
end

%% more problem dimensions
nx    = nb+nv + ng+nq;  %% number of standard OPF control variables
if nlin
  nz = size(mpc.A, 2) - nx; %% number of user z variables
  if nz < 0
    error('opf_setup: user supplied A matrix must have at least %d columns.', nx);
  end
else
  nz = 0;               %% number of user z variables
  if nw                 %% still need to check number of columns of N
    if size(mpc.N, 2) ~= nx;
      error('opf_setup: user supplied N matrix must have %d columns.', nx);
    end
  end
end

%% construct OPF model object
om = opf_model(mpc);
if ~isempty(pwl1)
  om.userdata.pwl1 = pwl1;
end
if dc
  %% user data
  om.userdata.Bf = Bf;
  om.userdata.Pfinj = Pfinj;
  om.userdata.iang = iang;

  %% optimization variables
  om.add_var('Va', nb, Va, Val, Vau);
  om.add_var('Pg', ng, Pg, Pmin, Pmax);

  %% linear constraints
  om.add_lin_constraint('Pmis', Amis, bmis, bmis, {'Va', 'Pg'});    %% nb
  om.add_lin_constraint('Pf',  Bf(il,:), -upt, upf, {'Va'});        %% nl2
  om.add_lin_constraint('ang', Aang, lang, uang, {'Va'});           %% nang

  %% quadratic generator costs
  if ~isempty(cpg)
    om.add_quad_cost('polPg', Qpg, cpg, kpg, {'Pg'});
  end
else
  %% user data
  om.userdata.Apqdata = Apqdata;
  om.userdata.iang = iang;

  %% optimization variables
  if vcart
      Vclim = 1.1 * bus(:, VMAX);
      om.add_var('Vr', nb, Vr, -Vclim, Vclim);
      om.add_var('Vi', nb, Vi, -Vclim, Vclim);
  else
      om.add_var('Va', nb, Va, Val, Vau);
      om.add_var('Vm', nb, Vm, bus(:, VMIN), bus(:, VMAX));
  end
  om.add_var('Pg', ng, Pg, Pmin, Pmax);
  om.add_var('Qg', ng, Qg, Qmin, Qmax);

  %% nonlinear constraints
  om.add_nln_constraint(mis_cons, [nb;nb], 1, fcn_mis, hess_mis, nodal_balance_vars);
  if legacy_formulation
    om.add_nln_constraint({'Sf', 'St'}, [nl;nl], 0, fcn_flow, hess_flow, flow_lim_vars);
  else
    om.add_nln_constraint({'Sf', 'St'}, [nl2;nl2], 0, fcn_flow, hess_flow, flow_lim_vars);
  end
  if vcart
    om.userdata.veq = veq;  %% buses with voltage magnitude equality constraints
    om.userdata.viq = viq;  %% buses with voltage magnitude limits
    om.add_nln_constraint('Vref', length(refs), 1, fcn_vref, hess_vref, {'Vr', 'Vi'});
    if nveq
      om.add_nln_constraint('Veq', nveq, 1, fcn_veq, hess_veq, {'Vr', 'Vi'});
    end
    om.add_nln_constraint({'Vmin', 'Vmax'}, [nvlims;nvlims], 0, fcn_vlim, hess_vlim, {'Vr', 'Vi'});
    om.add_nln_constraint({'angL', 'angU'}, [nang;nang], 0, fcn_ang, hess_ang, {'Vr', 'Vi'});
  end

  %% linear constraints
  om.add_lin_constraint('PQh', Apqh, [], ubpqh, {'Pg', 'Qg'});      %% npqh
  om.add_lin_constraint('PQl', Apql, [], ubpql, {'Pg', 'Qg'});      %% npql
  om.add_lin_constraint('vl',  Avl, lvl, uvl,   {'Pg', 'Qg'});      %% nvl
  if ~vcart
    om.add_lin_constraint('ang', Aang, lang, uang, {'Va'});         %% nang
  end

  %% polynomial generator costs
  if ~legacy_formulation
    %% quadratic/linear generator costs
    if ~isempty(cpg)
      om.add_quad_cost('polPg', Qpg, cpg, kpg, {'Pg'});
    end
    if ~isempty(cqg)
      om.add_quad_cost('polQg', Qqg, cqg, kqg, {'Qg'});
    end

    %% higher order polynomial generator costs
    if ~isempty(ip3)
      om.add_nln_cost('polPg', 1, cost_Pg, {'Pg'});
    end
    if ~isempty(qcost) && ~isempty(iq3)
      om.add_nln_cost('polQg', 1, cost_Qg, {'Qg'});
    end
  end
end

%% y vars, constraints for piece-wise linear gen costs
if ny > 0
  om.add_var('y', ny);
  om.add_lin_constraint('ycon', Ay, [], by, ycon_vars);             %% ncony
  if dc || ~legacy_formulation
    om.add_quad_cost('pwl', [], ones(ny, 1), 0, {'y'});
  end
end

%% add user vars, constraints and costs (as specified via A, ..., N, ...)
if nz > 0
  om.add_var('z', nz, z0, zl, zu);
  user_vars{end+1} = 'z';
end
if nlin
  om.add_lin_constraint('usr', mpc.A, lbu, ubu, user_vars);         %% nlin
end
if nnle
  for k = 1:length(mpc.user_constraints.nle)
    nlc = mpc.user_constraints.nle{k};
    fcn  = eval(['@(x)' nlc{3} '(x, nlc{6}{:})']);
    hess = eval(['@(x, lam)' nlc{4} '(x, lam, nlc{6}{:})']);
    om.add_nln_constraint(nlc{1:2}, 1, fcn, hess, nlc{5});
  end
end
if nnli
  for k = 1:length(mpc.user_constraints.nli)
    nlc = mpc.user_constraints.nli{k};
    fcn  = eval(['@(x)' nlc{3} '(x, nlc{6}{:})']);
    hess = eval(['@(x, lam)' nlc{4} '(x, lam, nlc{6}{:})']);
    om.add_nln_constraint(nlc{1:2}, 0, fcn, hess, nlc{5});
  end
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
  om.add_legacy_cost('usr', user_cost, user_vars);
end

%% execute userfcn callbacks for 'formulation' stage
om = run_userfcn(userfcn, 'formulation', om, mpopt);

%% implement legacy user costs using quadratic or general non-linear costs
cp = om.params_legacy_cost();   %% construct/fetch the parameters
[N, H, Cw, rh, mm] = deal(cp.N, cp.H, cp.Cw, cp.rh, cp.mm);
[nw, nx] = size(N);
if nw
    if any(cp.dd ~= 1) || any(cp.kk)    %% not simple quadratic form
        if dc                           %% (includes "dead zone" or
            if any(cp.dd ~= 1)          %%  quadratic "penalty")
                error('opf_setup: DC OPF can only handle legacy user-defined costs with d = 1');
            end
            if any(cp.kk)
                error('opf_setup: DC OPF can only handle legacy user-defined costs with no "dead zone", i.e. k = 0');
            end
        elseif ~legacy_formulation
            %% use general nonlinear cost to implement legacy user cost
            legacy_cost_fcn = @(x)opf_legacy_user_cost_fcn(x, cp);
            om.add_nln_cost('usr', 1, legacy_cost_fcn);
        end
    else                                %% simple quadratic form
        %% use a quadratic cost to implement legacy user cost
        if dc || ~legacy_formulation
            %% f = 1/2 * w'*H*w + Cw'*w, where w = diag(mm)*(N*x - rh)
            %% Let: MN = diag(mm)*N
            %%      MR = M * rh
            %%      HMR  = H  * MR;
            %%      HtMR = H' * MR;
            %%  =>   w = MN*x - MR
            %% f = 1/2 * (MN*x - MR)'*H*(MN*x - MR) + Cw'*(MN*x - MR)
            %%   = 1/2 * x'*MN'*H*MN*x +
            %%          (Cw'*MN - 1/2 * MR'*(H+H')*MN)*x +
            %%          1/2 * MR'*H*MR - Cw'*MR
            %%   = 1/2 * x'*Q*w + c'*x + k
    
            [N, H, Cw, rh, mm] = deal(cp.N, cp.H, cp.Cw, cp.rh, cp.mm);
            nw = size(N, 1);            %% number of general cost vars, w
            M    = sparse(1:nw, 1:nw, mm, nw, nw);
            MN   = M * N;
            MR   = M * rh;
            HMR  = H  * MR;
            HtMR = H' * MR;
            Q = MN' * H * MN;
            c = full(MN' * (Cw - 1/2*(HMR+HtMR)));
            k = (1/2 * HtMR - Cw)' * MR;
            om.add_quad_cost('usr', Q, c, k);
        end
    end
end

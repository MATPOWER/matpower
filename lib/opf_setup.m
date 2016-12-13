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

if dc
  %% ignore reactive costs for DC
  mpc.gencost = pqcost(mpc.gencost, ng);

  %% reduce A and/or N from AC dimensions to DC dimensions, if needed
  if nusr || nw
    acc = [nb+(1:nb) 2*nb+ng+(1:ng)];   %% Vm and Qg columns
    if nusr && size(mpc.A, 2) >= 2*nb + 2*ng
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
    elseif use_vg > 0 && use_vg < 1     %% fractional value
        %% use weighted avg between original Vmin/Vmax limits and Vg
        mpc.bus(ib, VMAX) = (1-use_vg) * mpc.bus(ib, VMAX) + use_vg * Vmax(ib);
        mpc.bus(ib, VMIN) = (1-use_vg) * mpc.bus(ib, VMIN) + use_vg * Vmin(ib);
    else
        error('opf_setup: option ''opf.use_vg'' (= %g) cannot be negative or greater than 1', use_vg);
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
Vm   = bus(:, VM);
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

  %% dispatchable load, constant power factor constraints
  [Avl, lvl, uvl]  = makeAvl(baseMVA, gen);
  
  %% generator PQ capability curve constraints
  [Apqh, ubpqh, Apql, ubpql, Apqdata] = makeApq(baseMVA, gen);

  user_vars = {'Va', 'Vm', 'Pg', 'Qg'};
  ycon_vars = {'Pg', 'Qg', 'y'};
end

%% voltage angle reference constraints
Vau = Inf(nb, 1);
Val = -Vau;
Vau(refs) = Va(refs);
Val(refs) = Va(refs);

%% branch voltage angle difference limits
[Aang, lang, uang, iang]  = makeAang(baseMVA, branch, nb, mpopt);

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
if nusr
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
  om = userdata(om, 'pwl1', pwl1);
end
if dc
  om = userdata(om, 'Bf', Bf);
  om = userdata(om, 'Pfinj', Pfinj);
  om = userdata(om, 'iang', iang);
  om = add_vars(om, 'Va', nb, Va, Val, Vau);
  om = add_vars(om, 'Pg', ng, Pg, Pmin, Pmax);
  om = add_constraints(om, 'Pmis', Amis, bmis, bmis, {'Va', 'Pg'}); %% nb
  om = add_constraints(om, 'Pf',  Bf(il,:), -upt, upf, {'Va'});     %% nl2
  om = add_constraints(om, 'ang', Aang, lang, uang, {'Va'});        %% nang
else
  om = userdata(om, 'Apqdata', Apqdata);
  om = userdata(om, 'iang', iang);
  om = add_vars(om, 'Va', nb, Va, Val, Vau);
  om = add_vars(om, 'Vm', nb, Vm, bus(:, VMIN), bus(:, VMAX));
  om = add_vars(om, 'Pg', ng, Pg, Pmin, Pmax);
  om = add_vars(om, 'Qg', ng, Qg, Qmin, Qmax);
  om = add_constraints(om, 'Pmis', nb, 'nonlinear');
  om = add_constraints(om, 'Qmis', nb, 'nonlinear');
  om = add_constraints(om, 'Sf', nl, 'nonlinear');
  om = add_constraints(om, 'St', nl, 'nonlinear');
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

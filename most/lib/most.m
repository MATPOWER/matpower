function mdo = most(mdi, mpopt)
%MOST MATPOWER Optimal Scheduling Tool
%   MDO = MOST(MDI)
%   MDO = MOST(MDI, MPOPT)
%
%   Solves a multiperiod, stochastic, contingency constrained, optimal
%   power flow problem with linear constraints and unit commitment.
%   Depending on inputs it may include DC power flow constraints or
%   a simple total power balance condition.
%
%   Inputs:
%       MDI   MOST data structure, input
%           (see MOST User's Manual for details)
%       MPOPT   MATPOWER options struct, relevant fields are (default
%               value in parens):
%           verbose - see 'help mpoption'
%           <solver specific options> - e.g. cplex, gurobi, etc,
%                     see 'help mpoption'
%           most.build_model (1) - build the MIQP, both constraints and
%                   standard costs (not coordination cost) and store in
%                   QP field of MDO
%           most.solve_model (1) - solve the MIQP; if coordination
%                   cost exists, update it; requires either 'most.build_model'
%                   set to 1 or MDI.QP must contain previously built model
%           most.resolve_new_cost (0) - use when MIQP is already built and
%                   unchanged except for new coordination cost
%           most.dc_model (1) - use DC flow network model as opposed to simple
%                   generation = demand constraint
%           most.fixed_res (-1) - include fixed zonal reserve contstraints,
%                   -1 = if present, 1 = always include, 0 = never include
%           most.q_coordination (0) - create Qg variables for reactive power
%                   coordination
%           most.security_constraints (-1) - include contingency contstraints,
%                   -1 = if present, 1 = always include, 0 = never include
%           most.storage.terminal_target (-1) - constrain the expected terminal
%                   storage to target value, if present (1 = always, 0 = never)
%           most.storage.cyclic (0) - if 1, then initial storage is a variable
%                   constrained to = final expected storage; can't be
%                   simultaneously true with most.storage.terminal_target
%           most.uc.run (-1) - flag to indicate whether to perform unit
%                   commitment; 0 = do NOT perform UC, 1 = DO perform UC,
%                   -1 = perform UC if MDI.UC.CommitKey is present/non-empty
%           most.uc.cyclic (0) - commitment restrictions (e.g. min up/down
%                   times) roll over from end of horizon back to beginning
%           most.alpha (0) - 0 = contingencies happen at beginning of period,
%                   1 = at end of period
%           most.solver ('DEFAULT') - see ALG argument to MIQPS_MATPOWER or
%                   QPS_MATPOWER for details
%           most.skip_prices (0) - skip price computation stage for mixed
%                   integer problems, see 'help miqps_matpower' for details
%           most.price_stage_warn_tol (1e-7) - tolerance on the objective fcn
%               value and primal variable relative match required to avoid
%               mis-match warning message, see 'help miqps_matpower' for details
%
%   Outputs:
%       MDO   MOST data structure, output
%           (see MOST User's Manual for details)


%   MOST
%   Copyright (c) 2010-2018, Power Systems Engineering Research Center (PSERC)
%   by Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Nacional de Colombia
%   and Ray Zimmerman, PSERC Cornell
%
%   This file is part of MOST.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/most for more info.

t0 = tic;

%% default arguments
if nargin < 2
    mpopt = mpoption;       %% use default options
end

verbose = mpopt.verbose;

if verbose
    fprintf('\n=============================================================================\n');
    fprintf(  '          MATPOWER Optimal Scheduling Tool  --  MOST Version %s\n', mostver());
    fprintf(  '          A multiperiod stochastic secure OPF with unit commitment\n');
    fprintf(  '                       -----  Built on MATPOWER  -----\n');
    fprintf(  '  by Carlos E. Murillo-Sanchez, Universidad Nacional de Colombia--Manizales\n');
    fprintf(  '                  and Ray D. Zimmerman, Cornell University\n');
    fprintf(  '       (c) 2012-2018 Power Systems Engineering Research Center (PSERC)       \n');
    fprintf(  '=============================================================================\n');
end

%% if you want to do a normal solve, you have to create the QP
if mpopt.most.solve_model && ~mpopt.most.resolve_new_cost
  mpopt = mpoption(mpopt, 'most.build_model', 1);
end
if ~mpopt.most.build_model && ~mpopt.most.solve_model
  error('most: Ah ... are you sure you want to do nothing? (either ''most.build_model'' or ''most.solve_model'' must be true)');
end

%% set up some variables we use throughout
ng = size(mdi.mpc.gen, 1);
nt = mdi.idx.nt;
ns = length(mdi.Storage.UnitIdx);
mdi.idx.ng = ng;
mdi.idx.ns = ns;
baseMVA = mdi.mpc.baseMVA;
for t = 1:nt
  mdi.idx.nj(t) = length(mdi.tstep(t).OpCondSched);
end
if ~isfield(mdi, 'UC') || ~isfield(mdi.UC, 'CommitSched') || ...
        isempty(mdi.UC.CommitSched)
  if isfield(mdi, 'CommitSched') && ~isempty(mdi.CommitSched)
    warning('-----  most: MDI.CommitSched has moved to MDI.UC.CommitSched, please update your code.  -----');
    mdi.UC.CommitSched = mdi.CommitSched;
  else
    error('most: commitment schedule must be provided in MSPD_IN.UC.CommitSched');
  end
end

% set up model options
UC = mpopt.most.uc.run;
if UC == -1
  if isempty(mdi.UC.CommitKey)
    UC = 0;
  else
    UC = 1;
  end
end
mo = struct(...
    'DCMODEL',                      mpopt.most.dc_model, ...
    'IncludeFixedReserves',         mpopt.most.fixed_res, ...
    'SecurityConstrained',          mpopt.most.security_constraints, ...
    'QCoordination',                mpopt.most.q_coordination, ...
    'ForceCyclicStorage',           mpopt.most.storage.cyclic, ...
    'CyclicCommitment',             mpopt.most.uc.cyclic, ...
    'ForceExpectedTerminalStorage', mpopt.most.storage.terminal_target, ...
    'alpha',                        mpopt.most.alpha ...
);
if mo.IncludeFixedReserves == -1
  if isfield(mdi, 'FixedReserves') && isfield(mdi.FixedReserves(1,1,1), 'req')
    mo.IncludeFixedReserves = 1;
  else
    mo.IncludeFixedReserves = 0;
  end
end
if mo.SecurityConstrained == -1
  if isfield(mdi, 'cont') && isfield(mdi.cont(1,1), 'contab') && ...
          ~isempty(mdi.cont(1,1).contab)
    mo.SecurityConstrained = 1;
  else
    mo.SecurityConstrained = 0;
  end
end
if mo.ForceExpectedTerminalStorage == -1
  if isempty(mdi.Storage.ExpectedTerminalStorageAim) && ...
     isempty(mdi.Storage.ExpectedTerminalStorageMin) && ...
     isempty(mdi.Storage.ExpectedTerminalStorageMax)
    mo.ForceExpectedTerminalStorage = 0;
  else
    mo.ForceExpectedTerminalStorage = 1;
  end
end
if mo.IncludeFixedReserves && ~(isfield(mdi, 'FixedReserves') && ...
      isfield(mdi.FixedReserves(1,1,1), 'req'))
  error('most: MDI.FixedReserves(t,j,k) must be specified when MPOPT.most.fixed_res = 1');
end
if mo.SecurityConstrained && ~(isfield(mdi, 'cont') && ...
      isfield(mdi.cont(1,1), 'contab') && ~isempty(mdi.cont(1,1).contab))
  error('most: MDI.cont(t,j).contab cannot be empty when MPOPT.most.security_constraints = 1');
end
if mo.IncludeFixedReserves && mo.SecurityConstrained
  warning('most: Using MPOPT.most.fixed_res = 1 and MPOPT.most.security_constraints = 1 together is not recommended.');
end
if mo.ForceExpectedTerminalStorage == 1;
  if mo.ForceCyclicStorage
    error('most: storage model cannot be both cyclic and include a terminal target value; must change MPOPT.most.storage.cyclic or MPOPT.most.storage.terminal_target');
  end
  if ns && isempty(mdi.Storage.ExpectedTerminalStorageAim) && ...
           isempty(mdi.Storage.ExpectedTerminalStorageMin) && ...
           isempty(mdi.Storage.ExpectedTerminalStorageMax)
    error('most: MDI.Storage.ExpectedTerminalStorageAim|Min|Max cannot all be empty when MPOPT.most.storage.terminal_target = 1');
  end
  if ~isempty(mdi.Storage.ExpectedTerminalStorageAim)
    mdi.Storage.ExpectedTerminalStorageMin = mdi.Storage.ExpectedTerminalStorageAim;
    mdi.Storage.ExpectedTerminalStorageMax = mdi.Storage.ExpectedTerminalStorageAim;
  end
end
if UC && (~isfield(mdi.UC, 'CommitKey') || isempty(mdi.UC.CommitKey))
  error('most: cannot run unit commitment without specifying MDI.UC.CommitKey');
end

if ns
  if isempty(mdi.Storage.InitialStorage)
    error('most: Storage.InitialStorage must be specified');
  end
  if isempty(mdi.Storage.TerminalChargingPrice0)
    mdi.Storage.TerminalChargingPrice0 = mdi.Storage.TerminalStoragePrice;
  end
  if isempty(mdi.Storage.TerminalDischargingPrice0)
    mdi.Storage.TerminalDischargingPrice0 = mdi.Storage.TerminalStoragePrice;
  end
  if isempty(mdi.Storage.TerminalChargingPriceK)
    mdi.Storage.TerminalChargingPriceK = mdi.Storage.TerminalStoragePrice;
  end
  if isempty(mdi.Storage.TerminalDischargingPriceK)
    mdi.Storage.TerminalDischargingPriceK = mdi.Storage.TerminalStoragePrice;
  end
  if isempty(mdi.Storage.MinStorageLevel)
    error('most: Storage.MinStorageLevel must be specified');
  else
    MinStorageLevel = mdi.Storage.MinStorageLevel;
  end
  if size(MinStorageLevel, 1) == 1 && ns > 1    %% expand rows
    MinStorageLevel = ones(ns, 1) * MinStorageLevel;
  end
  if size(MinStorageLevel, 2) == 1 && nt > 1    %% expand cols
    MinStorageLevel = MinStorageLevel * ones(1, nt);
  end
  if isempty(mdi.Storage.MaxStorageLevel)
    error('most: Storage.MaxStorageLevel must be specified');
  else
    MaxStorageLevel = mdi.Storage.MaxStorageLevel;
  end
  if size(MaxStorageLevel, 1) == 1 && ns > 1    %% expand rows
    MaxStorageLevel = ones(ns, 1) * MaxStorageLevel;
  end
  if size(MaxStorageLevel, 2) == 1 && nt > 1    %% expand cols
    MaxStorageLevel = MaxStorageLevel * ones(1, nt);
  end
  if isempty(mdi.Storage.InEff)
    InEff = 1;                      %% no efficiency loss by default
  else
    InEff = mdi.Storage.InEff;
  end
  if size(InEff, 1) == 1 && ns > 1  %% expand rows
    InEff = ones(ns, 1) * InEff;
  end
  if size(InEff, 2) == 1 && nt > 1  %% expand cols
    InEff = InEff * ones(1, nt);
  end
  if isempty(mdi.Storage.OutEff)
    OutEff = 1;                     %% no efficiency loss by default
  else
    OutEff = mdi.Storage.OutEff;
  end
  if size(OutEff, 1) == 1 && ns > 1 %% expand rows
    OutEff = ones(ns, 1) * OutEff;
  end
  if size(OutEff, 2) == 1 && nt > 1 %% expand cols
    OutEff = OutEff * ones(1, nt);
  end
  if isempty(mdi.Storage.LossFactor)
    LossFactor = 0;                     %% no losses by default
  else
    LossFactor = mdi.Storage.LossFactor;
  end
  if size(LossFactor, 1) == 1 && ns > 1 %% expand rows
    LossFactor = ones(ns, 1) * LossFactor;
  end
  if size(LossFactor, 2) == 1 && nt > 1 %% expand cols
    LossFactor = LossFactor * ones(1, nt);
  end
  if isempty(mdi.Storage.rho)
    rho = 1;                        %% use worst case by default (for backward compatibility)
  else
    rho = mdi.Storage.rho;
  end
  if size(rho, 1) == 1 && ns > 1    %% expand rows
    rho = ones(ns, 1) * rho;
  end
  if size(rho, 2) == 1 && nt > 1    %% expand cols
    rho = rho * ones(1, nt);
  end
  if isempty(mdi.Storage.InitialStorageLowerBound)
    if mo.ForceCyclicStorage        %% lower bound for var s0, take from t=1
      mdi.Storage.InitialStorageLowerBound = MinStorageLevel(:, 1);
    else                            %% Sm(0), default = fixed param s0
      mdi.Storage.InitialStorageLowerBound = mdi.Storage.InitialStorage;
    end
  end
  if isempty(mdi.Storage.InitialStorageUpperBound)
    if mo.ForceCyclicStorage        %% upper bound for var s0, take from t=1
      mdi.Storage.InitialStorageUpperBound = MaxStorageLevel(:, 1);
    else                            %% Sp(0), default = fixed param s0
      mdi.Storage.InitialStorageUpperBound = mdi.Storage.InitialStorage;
    end
  end
  
  LossCoeff = mdi.Delta_T * LossFactor/2;
  beta1 = (1-LossCoeff) ./ (1+LossCoeff);
  beta2 = 1 ./ (1+LossCoeff);
  beta3 = 1 ./ (1/(1-mo.alpha) + LossCoeff);
  beta4 = mo.alpha/(1-mo.alpha) * beta2 .* beta3;
  beta5 = beta1 ./ beta2 .* (beta3 + beta4);
  beta2EtaIn      = mdi.Delta_T * beta2 .* InEff;
  beta2overEtaOut = mdi.Delta_T * beta2 ./ OutEff;
  beta3EtaIn      = mdi.Delta_T * beta3 .* InEff;
  beta3overEtaOut = mdi.Delta_T * beta3 ./ OutEff;
  beta4EtaIn      = mdi.Delta_T * beta4 .* InEff;
  beta4overEtaOut = mdi.Delta_T * beta4 ./ OutEff;
  diagBeta2EtaIn1      = spdiags(beta2EtaIn(:,1),      0, ns, ns);
  diagBeta2overEtaOut1 = spdiags(beta2overEtaOut(:,1), 0, ns, ns);
  diagBeta3EtaIn1      = spdiags(beta3EtaIn(:,1),      0, ns, ns);
  diagBeta3overEtaOut1 = spdiags(beta3overEtaOut(:,1), 0, ns, ns);
  diagBeta4EtaIn1      = spdiags(beta4EtaIn(:,1),      0, ns, ns);
  diagBeta4overEtaOut1 = spdiags(beta4overEtaOut(:,1), 0, ns, ns);
end
if ~isfield(mdi.idx, 'ntds') || isempty(mdi.idx.ntds) || ~mdi.idx.ntds
  ntds = 0;
  nzds = 0;
  nyds = 0;
else
  ntds = mdi.idx.ntds;
  nzds = size(mdi.dstep(1).A, 1);
  nyds = size(mdi.dstep(1).D, 1);   % # of outputs of dynamical system
end
mdi.idx.ntds = ntds;
mdi.idx.nzds = nzds;
mdi.idx.nyds = nyds;

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
[CT_LABEL, CT_PROB, CT_TABLE, CT_TBUS, CT_TGEN, CT_TBRCH, CT_TAREABUS, ...
    CT_TAREAGEN, CT_TAREABRCH, CT_ROW, CT_COL, CT_CHGTYPE, CT_REP, ...
    CT_REL, CT_ADD, CT_NEWVAL, CT_TLOAD, CT_TAREALOAD, CT_LOAD_ALL_PQ, ...
    CT_LOAD_FIX_PQ, CT_LOAD_DIS_PQ, CT_LOAD_ALL_P, CT_LOAD_FIX_P, ...
    CT_LOAD_DIS_P, CT_TGENCOST, CT_TAREAGENCOST, CT_MODCOST_F, ...
    CT_MODCOST_X] = idx_ct;

%% check that bus numbers are equal to indices to bus (one set of bus numbers)
nb  = size(mdi.mpc.bus, 1);
if any(mdi.mpc.bus(:, BUS_I) ~= (1:nb)')
    error('most: buses must be numbered consecutively in bus matrix; use ext2int() to convert to internal ordering')
end

% Make data tables with full # of cols and add also pseudo OPF results to
% be able to run printpf on them
mdi.mpc.bus(:, MU_VMIN) = 0;
mdi.mpc.gen(:, MU_QMIN) = 0;
mdi.mpc.branch(:,MU_ANGMAX) = 0;
mdi.mpc.f = 0;
mdi.mpc.et = 0;
mdi.mpc.success = 1;

if mpopt.most.build_model
  if verbose
    fprintf('- Building indexing structures.\n');
  end

  %% save model options in data structure
  mdi.DCMODEL                               = mo.DCMODEL;
  mdi.IncludeFixedReserves                  = mo.IncludeFixedReserves;
  mdi.SecurityConstrained                   = mo.SecurityConstrained;
  mdi.QCoordination                         = mo.QCoordination;
  mdi.Storage.ForceCyclicStorage            = mo.ForceCyclicStorage;
  mdi.Storage.ForceExpectedTerminalStorage  = mo.ForceExpectedTerminalStorage;
  mdi.UC.run                                = UC;
  mdi.UC.CyclicCommitment                   = mo.CyclicCommitment;
  mdi.alpha                                 = mo.alpha;
  if ~isfield(mdi, 'OpenEnded'), mdi.OpenEnded = 1; end

  if UC
    % Make sure MinUp and MinDown are all >= 1
    if any(mdi.UC.MinUp < 1) && any(mdi.UC.MinUp < 1)
        error('most: UC.MinUp and UC.MinDown must all be >= 1');
    end
    % Unless something is forced off in mdi.CommitKey, or as a result of
    % not fulfilling its mdi.UC.MinDown in early periods, it should be available
    % for commitment and thus a contingency including its outage should not
    % be deleted.
    mdi.UC.CommitSched = (mdi.UC.CommitKey >= 0);   % Treat anything but -1 as on.
    if ~mdi.UC.CyclicCommitment
      for i = 1:ng
        if mdi.UC.InitialState(i) < 0
          nn = mdi.UC.MinDown(i) + mdi.UC.InitialState(i);  % time to go before startup
          if nn > 0
            mdi.UC.CommitSched(i, 1:nn) = 0;
          end
        elseif mdi.UC.InitialState(i) > 0
          nn = mdi.UC.MinUp(i) - mdi.UC.InitialState(i);    % time to go before shutdown
          if nn > 0
            mdi.UC.CommitSched(i, 1:nn) = 1;
          end
        end
      end
    end
  end
  % From now on, mdi.UC.CommitSched has zeros for stuff that is definitely
  % off, and ones for stuff that might be on, so those zeros can be used to
  % trim off contingencies that won't happen.
  % Start by creating the base flow data for all scenarios
  for t = 1:nt
    for j = 1:mdi.idx.nj(t)
      mpc = mdi.mpc;
      mpc.gen(:, GEN_STATUS) = mdi.UC.CommitSched(:, t);

      %% for backward compatibility, putting time dependent energy offer
      %% data in offer(t).gencost is deprecated, please use profiles
      if isfield(mdi.offer(t), 'gencost') && ~isempty(mdi.offer(t).gencost)
        mpc.gencost = mdi.offer(t).gencost;
      end

      if ~isempty(mdi.tstep(t).OpCondSched(j).tab)
        changelist = unique(mdi.tstep(t).OpCondSched(j).tab(:, CT_LABEL));
        for label = changelist'
          mpc = apply_changes(label, mpc, mdi.tstep(t).OpCondSched(j).tab);
        end
      end
      mdi.flow(t,j,1).mpc = mpc;
      mdi.idx.nb(t,j,1) = size(mdi.flow(t,j,1).mpc.bus, 1);
      mdi.idx.ny(t,j,1) = length(find(mdi.flow(t,j,1).mpc.gencost(:, MODEL) == PW_LINEAR));
    end
  end
  % Then continue to create contingent flow scenarios, deleting any
  % redundant contingencies (i.e., decommitting a gen or branch when its
  % status is guaranteed to be off). No rows are deleted from gen or branch,
  % but the number of contingencies can indeed change.
  for t = 1:nt
    % Set default ramp reserve mask, if not provided
    if ~isfield(mdi.tstep(t), 'TransMask') || isempty(mdi.tstep(t).TransMask)
      mdi.tstep(t).TransMask = ones(size(mdi.tstep(t).TransMat));
    end
    % First get current step's scenario probabilities
    if t == 1
      scenario_probs = mdi.tstep(1).TransMat; % the probability of the initial state is 1
    else
      scenario_probs = mdi.tstep(t).TransMat * mdi.CostWeights(1, 1:mdi.idx.nj(t-1), t-1)'; % otherwise compute from previous step base cases
    end
    mdi.StepProb(t) = sum(scenario_probs); % probability of making it to the t-th step
    if mdi.SecurityConstrained
      for j = 1:mdi.idx.nj(t)
        [tmp, ii] = sort(mdi.cont(t,j).contab(:, CT_LABEL)); %sort in ascending contingency label
        contab = mdi.cont(t,j).contab(ii, :);
        rowdecomlist = ones(size(contab,1), 1);
        for l = 1:size(contab, 1)
          if contab(l, CT_TABLE) == CT_TGEN  && contab(l, CT_COL) == GEN_STATUS ...
              && contab(l, CT_CHGTYPE) == CT_REP && contab(l, CT_NEWVAL) == 0 ... % gen turned off
              && mdi.flow(t,j,1).mpc.gen(contab(l, CT_ROW), GEN_STATUS) <= 0    % but it was off on input
           rowdecomlist(l) = 0;
          elseif contab(l, CT_TABLE) == CT_TBRCH && contab(l, CT_COL) == BR_STATUS ...
              && contab(l, CT_CHGTYPE) == CT_REP && contab(l, CT_NEWVAL) == 0 ... % branch taken out
              && mdi.flow(t,j,1).mpc.branch(contab(l, CT_ROW), BR_STATUS) <= 0  % but it was off on input
            rowdecomlist(l) = 0;
          end
        end
        contab = contab(rowdecomlist ~= 0, :);
        mdi.cont(t, j).contab = contab;
        clist = unique(contab(:, CT_LABEL));
        mdi.idx.nc(t, j) = length(clist);
        k = 2;
        for label = clist'
          mdi.flow(t, j, k).mpc = apply_changes(label, mdi.flow(t, j, 1).mpc, contab);
          ii = find( label == contab(:, CT_LABEL) );
          mdi.CostWeights(k, j, t) = contab(ii(1), CT_PROB);
          mdi.idx.nb(t, j, k) = size(mdi.flow(t, j, k).mpc.bus, 1);
          mdi.idx.ny(t, j, k) = length(find(mdi.flow(t, j, 1).mpc.gencost(:, MODEL) == PW_LINEAR));
          k = k + 1;
        end
        mdi.CostWeights(1, j, t) = 1 - sum(mdi.CostWeights(2:mdi.idx.nc(t,j)+1, j, t));
        mdi.CostWeights(1:mdi.idx.nc(t,j)+1, j, t) = scenario_probs(j) * mdi.CostWeights(1:mdi.idx.nc(t,j)+1, j, t);
      end
    else
      for j = 1:mdi.idx.nj(t)
        mdi.idx.nc(t, j) = 0;
        mdi.CostWeights(1, j, t) = scenario_probs(j);
      end
    end
  end

  % Compute adjusted (for alpha) cost weights for objective function
  if mdi.SecurityConstrained && mdi.alpha ~= 0
    for t = 1:nt
      for j = 1:mdi.idx.nj(t)
        mdi.CostWeightsAdj(1, j, t) = mdi.CostWeights(1, j, t);
        for k = 2:mdi.idx.nc(t,j)+1
          mdi.CostWeightsAdj(k, j, t) = (1-mdi.alpha) * mdi.CostWeights(k, j, t);
          mdi.CostWeightsAdj(1, j, t) = mdi.CostWeightsAdj(1, j, t) + mdi.alpha * mdi.CostWeights(k, j, t);
        end
      end
    end
  else
    mdi.CostWeightsAdj = mdi.CostWeights;
  end

  % If UC, also need to (possibly) modify gencosts so that each fm(p) at
  % p=0 is zero, so that fm(p) + u*c00 is equal to the original f(p). This
  % portion of the cost is valid at t if the unit is commited there, but
  % only for base scenarios and contingencies in which this particular unit
  % is not ousted, so must be careful later when probability-weighting the
  % corresponding u(i,t)!
  if UC
    if ~isfield(mdi.UC, 'c00') || isempty(mdi.UC.c00)   % if not empty assume
      mdi.UC.c00 = zeros(ng, nt);                       % contains correct info!
      for t = 1:nt
        for j = 1:mdi.idx.nj(t)
          for k = 1:mdi.idx.nc(t,j)+1
            mpc = mdi.flow(t,j,k).mpc;
            c00tjk = totcost(mpc.gencost, zeros(ng,1));
            mdi.UC.c00(:, t) = mdi.UC.c00(:, t) + mdi.CostWeightsAdj(k, j, t) * c00tjk;
            c0col = COST + mpc.gencost(:,NCOST) - 1;
            ipoly = find(mpc.gencost(:, MODEL) == POLYNOMIAL);
            ipwl  = find(mpc.gencost(:, MODEL) == PW_LINEAR);
            ii = sub2ind(size(mpc.gencost), ipoly, c0col(ipoly));
            mpc.gencost(ii) = mpc.gencost(ii) - c00tjk(ipoly);
            for i = ipwl'
              jj = COST+1:2:COST+2*mpc.gencost(i,NCOST)-1;
              mpc.gencost(i, jj) = mpc.gencost(i, jj) - c00tjk(i);
            end
            mpc.fixed_gencost = c00tjk;
            mdi.flow(t,j,k).mpc = mpc;
          end
        end
      end
    end
  end
  
  % Build variable indexing mechanism
  % Find total number of flows, buses and ny variables; (including offline gens)
  mdi.idx.nf_total = 0;
  mdi.idx.nb_total = 0;
  for t = 1:nt
    for j = 1:mdi.idx.nj(t)
      mdi.idx.nf_total = mdi.idx.nf_total + (1 + mdi.idx.nc(t,j));
      for k = 1:mdi.idx.nc(t,j)+1
        mdi.idx.nb_total = mdi.idx.nb_total + size(mdi.flow(t, j, k).mpc.bus, 1);
        ii = find(mdi.flow(t,j,k).mpc.gencost(:, MODEL) == PW_LINEAR);
        mdi.idx.ny(t,j,k) = length(ii);
      end
    end
  end
  mdi.idx.ns_total = ns * mdi.idx.nf_total;
  % Variable order resembles that of several C3SOPFs stacked together,
  % including the internally generated y variables, and then all of the
  % other new variables that are specific to HP, but excluding qg on one hand, and pc,
  % rp and rm since these now are common across several scenarios. So create first a matrix
  % of indices to the beginning of each c3sopf cell's vars.  Include the
  % mechanism for adding theta variables if we want to create DC flow restrictions.
  % Then start assigning the start and end indices for variables in each
  % c3sopf cell
  om = opt_model;
  nj_max = max(mdi.idx.nj);
  nc_max = max(max(mdi.idx.nc));
  Ing = speye(ng);
  Ins = speye(ns);
  if mdi.DCMODEL
    om.init_indexed_name('var', 'Va', {nt, nj_max, nc_max+1});
  end
  om.init_indexed_name('var', 'Pg', {nt, nj_max, nc_max+1});
  om.init_indexed_name('var', 'dPp', {nt, nj_max, nc_max+1});
  om.init_indexed_name('var', 'dPm', {nt, nj_max, nc_max+1});
  om.init_indexed_name('var', 'y', {nt, nj_max, nc_max+1});
  if mdi.IncludeFixedReserves
    om.init_indexed_name('var', 'R', {nt, nj_max, nc_max+1});
  end
  for t = 1:nt
    for j = 1:mdi.idx.nj(t)
      % first all angles if using DCMODEL
      for k = 1:mdi.idx.nc(t,j)+1
        if mdi.DCMODEL
          iref = find(mdi.flow(t,j,k).mpc.bus(:,BUS_TYPE) == REF);
          if verbose && length(iref) > 1
            errstr = ['\nmost: Warning: Multiple reference buses.\n', ...
              '           For a system with islands, a reference bus in each island\n', ...
              '           may help convergence, but in a fully connected system such\n', ...
              '           a situation is probably not reasonable.\n\n' ];
            fprintf(errstr);
          end
          Va0 = mdi.flow(t,j,k).mpc.bus(:,VA)*pi/180;
          Va_max = Inf(mdi.idx.nb(t,j,k), 1);
          Va_min = -Va_max;
          Va_min(iref) = mdi.flow(t,j,k).mpc.bus(iref,VA)*pi/180;
          Va_max(iref) = Va_min(iref);
          
          om.add_var('Va', {t,j,k}, mdi.idx.nb(t,j,k), Va0, Va_min, Va_max);
        end
      end
      % All active injections in c3sopf cell
      for k = 1:mdi.idx.nc(t,j)+1
        mpc = mdi.flow(t,j,k).mpc;
        genmask = mpc.gen(:,GEN_STATUS) > 0;
        p0 = genmask .* mpc.gen(:,PG) / baseMVA;
        if UC       % relax bounds here, enforced by uPmax, uPmin constraints
          pmin = genmask .* (min(mpc.gen(:, PMIN) / baseMVA, 0) - 1);
          pmax = genmask .* (max(mpc.gen(:, PMAX) / baseMVA, 0) + 1);
        else        % enforce bounds here, subject to flow's GEN_STATUS
          pmin = genmask .* mpc.gen(:, PMIN) / baseMVA;
          pmax = genmask .* mpc.gen(:, PMAX) / baseMVA;
        end
        om.add_var('Pg', {t,j,k}, ng, p0, pmin, pmax);
      end
      if mdi.IncludeFixedReserves
        for k = 1:mdi.idx.nc(t,j)+1
          mpc = mdi.flow(t,j,k).mpc;
          r = mdi.FixedReserves(t,j,k);
          nrz = size(r.req, 1); %% number of reserve zones
          if nrz > 1
              r.rgens = any(r.zones);   %% mask of gens available to provide reserves
          else
              r.rgens = r.zones;
          end
          r.igr = find(r.rgens);        %% indices of gens available to provide reserves
          ngr = length(r.igr);          %% number of gens available to provide reserves
          %% check data for consistent dimensions
          if size(r.zones, 1) ~= nrz
              error('most: the number of rows in FixedReserves(%d,%d,%d).req (%d) and FixedReserves(%d,%d,%d).zones (%d) must match', t, j, k, nrz, t, j, k, size(r.zones, 1));
          end
          if size(r.cost, 1) ~= ng && size(r.cost, 1) ~= ngr
              error('most: the number of rows in FixedReserves(%d,%d,%d).cost (%d) must equal the total number of generators (%d) or the number of generators able to provide reserves (%d)', t, j, k, size(r.cost, 1), ng, ngr);
          end
          if isfield(r, 'qty') && size(r.qty, 1) ~= size(r.cost, 1)
              error('most: FixedReserves(%d,%d,%d).cost (%d x 1) and FixedReserves(%d,%d,%d).qty (%d x 1) must be the same dimension', t, j, k, size(r.cost, 1), t, j, k, size(r.qty, 1));
          end
          %% convert both cost and qty from ngr x 1 to full ng x 1 vectors if necessary
          if size(r.cost, 1) < ng
              r.original.cost = r.cost;     %% save original
              cost = zeros(ng, 1);
              cost(r.igr) = r.cost;
              r.cost = cost;
              if isfield(r, 'qty')
                  r.original.qty = r.qty;   %% save original
                  qty = zeros(ng, 1);
                  qty(r.igr) = r.qty;
                  r.qty = qty;
              end
          end
          mdi.FixedReserves(t,j,k).rgens = r.rgens;
          mdi.FixedReserves(t,j,k).igr   = r.igr;
          if isfield(r, 'original')
              mdi.FixedReserves(t,j,k).original = r.original;
          end
          mdi.FixedReserves(t,j,k)       = r;   %% for cost & qty (now that fields match)
          Rmax = Inf(ngr, 1);               %% bound above by ...
          kk = find(mpc.gen(r.igr, RAMP_10));
          Rmax(kk) = mpc.gen(r.igr(kk), RAMP_10);   %% ... ramp rate and ...
          kk = find(r.qty(r.igr) < Rmax);
          Rmax(kk) = r.qty(r.igr(kk));      %% ... stated max reserve qty
          Rmax = Rmax / baseMVA;
          om.add_var('R', {t,j,k}, ngr, [], zeros(ngr, 1), Rmax);
        end
      end
      % All deltaP plus in c3sopf cell
      for k = 1:mdi.idx.nc(t,j)+1
        om.add_var('dPp', {t,j,k}, ng, [], zeros(ng,1), []);
      end
      % All deltaP minus in c3sopf cell
      for k = 1:mdi.idx.nc(t,j)+1
        om.add_var('dPm', {t,j,k}, ng, [], zeros(ng,1), []);
      end
      % All y variables in c3sopf cell - even if not committed.  There must
      % be a fixed cost associated with u(t,i,j) such that if u(t,i,j) = 0,
      % then the cost interpolated from the (x,y) pairs is zero, and if
      % u(t,i,j) = 1, then the fixed cost plus that interpolated from the
      % (x,y) pairs is as desired.
      for k = 1:mdi.idx.nc(t,j)+1
        om.add_var('y', {t,j,k}, mdi.idx.ny(t,j,k), [], [], []);
      end %
    end % for j
  end % for t
  % Continue with pc, rpp, rpm, one set for each time period
  om.init_indexed_name('var', 'Pc', {nt});
  om.init_indexed_name('var', 'Rpp', {nt});
  om.init_indexed_name('var', 'Rpm', {nt});
  for t = 1:nt
    om.add_var('Pc', {t}, ng);
    %% non-negativity on Rpp and Rpm is redundant, leave unbounded below
    %% (except where gen is off-line for all j and k)
    Rpmin = -Inf(ng,1);
    off = ones(ng,1);
    for j = 1:mdi.idx.nj(t);
      for k = 1:mdi.idx.nc(t,j)+1
        off = off & mdi.flow(t,j,k).mpc.gen(:, GEN_STATUS) <= 0;
      end
    end
    Rpmin(off == 1) = 0;
    om.add_var('Rpp', {t}, ng, [], Rpmin, mdi.offer(t).PositiveActiveReserveQuantity/baseMVA);
    om.add_var('Rpm', {t}, ng, [], Rpmin, mdi.offer(t).NegativeActiveReserveQuantity/baseMVA);
  end
  % Now load following ramping reserves.  In open ended problem, we need to
  % specify nt-1 ramping reserves, those needed to transition 1-2, 2-3, ..
  % (nt-1)-nt . The initial ramp constraint (from t=0 to t=1) is data, not a
  % variable.  But in terminal state (at t=nt+1) or cyclical problems (when
  % the t=nt to t=1 transition is also considered) we need nt ramping
  % reserves.
  if ~mdi.OpenEnded
    mdi.idx.ntramp = nt;
  else
    mdi.idx.ntramp = nt - 1;
  end
  om.init_indexed_name('var', 'Rrp', {mdi.idx.ntramp});
  om.init_indexed_name('var', 'Rrm', {mdi.idx.ntramp});
  for t = 1:mdi.idx.ntramp
    ramp30 = mdi.flow(t,1,1).mpc.gen(:,RAMP_30)*2*mdi.Delta_T;
    om.add_var('Rrp', {t}, ng, [], zeros(ng,1), ...
        min(mdi.offer(t).PositiveLoadFollowReserveQuantity, ramp30)/baseMVA);
  end
  for t = 1:mdi.idx.ntramp
    ramp30 = mdi.flow(t,1,1).mpc.gen(:,RAMP_30)*2*mdi.Delta_T;
    om.add_var('Rrm', {t}, ng, [], zeros(ng,1), ...
        min(mdi.offer(t).NegativeLoadFollowReserveQuantity, ramp30)/baseMVA);
  end
  % Continue with storage charge/discharge injections, one of each
  % for each flow; first all charge injections, then all discharge
  % injections
  om.init_indexed_name('var', 'Psc', {nt, nj_max, nc_max+1});
  om.init_indexed_name('var', 'Psd', {nt, nj_max, nc_max+1});
  for t = 1:nt
    for j = 1:mdi.idx.nj(t)
      for k = 1:mdi.idx.nc(t,j)+1
        if ns
          om.add_var('Psc', {t,j,k}, ns, [], [], zeros(ns,1));
        end
      end
    end
  end
  if ns
    for t = 1:nt
      for j = 1:mdi.idx.nj(t)
        for k = 1:mdi.idx.nc(t,j)+1
          om.add_var('Psd', {t,j,k}, ns, [], zeros(ns,1), []);
        end
      end
    end
  end
  % Continue with storage upper and lower bounds, one for each time period
  % and unit
  om.init_indexed_name('var', 'Sp', {nt});
  om.init_indexed_name('var', 'Sm', {nt});
  if ns
    for t = 1:nt
      om.add_var('Sp', {t}, ns, [], [], MaxStorageLevel(:,t)/baseMVA);
    end
    for t = 1:nt
      om.add_var('Sm', {t}, ns, [], MinStorageLevel(:,t)/baseMVA, []);
    end
  end
  % Possible initial storage quantities when using cyclic storage dispatch
  % so that initial storage = expected terminal storage is a constraint
  if ns && mdi.Storage.ForceCyclicStorage
    om.add_var('S0', ns, [], ...
        mdi.Storage.InitialStorageLowerBound / baseMVA, ...
        mdi.Storage.InitialStorageUpperBound / baseMVA);
  end
  % If there is a dynamical system with non-null state vector,
  % add those states here
  if nzds
    om.init_indexed_name('var', 'Z', {ntds});
    for t = 1:ntds
      if t == 1
        zmin = mdi.z1;
        zmax = mdi.z1;
      else
        zmin = mdi.dstep(t).zmin;
        zmax = mdi.dstep(t).zmax;
      end
      z0 = (zmax - zmin) / 2;
      om.add_var('Z', {t}, nzds, z0, zmin, zmax);
    end
  end
  % Now the integer variables; u variables mean on/off status
  if UC
    om.init_indexed_name('var', 'u', {nt});
    om.init_indexed_name('var', 'v', {nt});
    om.init_indexed_name('var', 'w', {nt});
    vt0 = char('B' * ones(1, ng));  % default variable type for u is binary
    for t = 1:nt
      umin = zeros(ng, 1);
      umax = ones(ng, 1);
      % min up/down restrictions on u
      if ~mdi.UC.CyclicCommitment
        % min up time has not passed yet since startup occured, force ON
        umin( (mdi.UC.InitialState > 0) & ...
              (t+mdi.UC.InitialState-mdi.UC.MinUp <= 0) ) = 1;
        % min down time has not passed yet since shutdown occured, force OFF
        umax( (mdi.UC.InitialState < 0) & ...
              (t-mdi.UC.InitialState-mdi.UC.MinDown <= 0) ) = 0;
      end
      % set limits for units forced ON or forced OFF
      iON = find(mdi.UC.CommitKey(:,t) == 2);
      iOFF = find(mdi.UC.CommitKey(:,t) == -1);
      umin(iON)  = 1;
      umax(iOFF) = 0;

      % set variable types
      vt = vt0;                 % initialize all variable types to binary
      vt(umin == umax) = 'C';   % make continuous for those that are fixed
      
      om.add_var('u', {t}, ng, zeros(ng, 1), umin, umax, vt);
    end
    % v variables mean startup events
    for t = 1:nt
      om.add_var('v', {t}, ng, zeros(ng, 1), zeros(ng, 1), ones(ng, 1));
    end
    % w variables mean shutdown events
    for t = 1:nt
      om.add_var('w', {t}, ng, zeros(ng, 1), zeros(ng, 1), ones(ng, 1));
    end
  end
  % An external program may be using coordination with AC flows, and in
  % that case we need corresponding Qg variables, whose only function is to
  % be constrained to zero if the commitment decision asks for that.
  if mdi.QCoordination
    om.init_indexed_name('var', 'Qg', {nt, nj_max, nc_max+1});
    for t = 1:nt
      for j = 1:mdi.idx.nj(t)
        for k = 1:mdi.idx.nc(t,j)+1
          mpc = mdi.flow(t,j,k).mpc;
          genmask = mpc.gen(:,GEN_STATUS) > 0;
          q0 = genmask .* mpc.gen(:, QG) / baseMVA;
          if UC     % relax bounds here, enforced by uQmax, uQmin constraints
            qmin = genmask .* (min(mpc.gen(:, QMIN) / baseMVA, 0) - 1);
            qmax = genmask .* (max(mpc.gen(:, QMAX) / baseMVA, 0) + 1);
          else      % enforce bounds here, subject to flow's GEN_STATUS
            qmin = genmask .* mpc.gen(:, QMIN) / baseMVA;
            qmax = genmask .* mpc.gen(:, QMAX) / baseMVA;
          end
          om.add_var('Qg', {t,j,k}, ng, q0, qmin, qmax);
        end
      end
    end
  end
  nvars = om.getN('var');
  mdi.idx.nvars = nvars;

  % Construct mechanism to keep track of expected storage states. This is
  % used for building some constraints in some types of problems, most
  % notably when there is a constraint on the terminal expected storage,
  % but it is also needed when there is a value associated with leftover
  % storage. Unfortunately this requires the construction of a substantial
  % mechanism for computing the expected terminal storage at any end-point
  % of the transition tree, be it a contingency or the terminal of the central
  % path at the end of the horizon.  Let SF(t) be the terminal storage value
  % at the end of the t-th period, assuming that we make it there without
  % contingencies, and let SI(t) be the initial amount of storage at that
  % same period. Then
  %     SI(t) = D(t) * SF(t-1).
  %     SF(t) = B1(t) * SI(t) + B2*[G(t)*x + H(t)]*x, and
  % Here, D(t) is created from the probability transition matrix, restricted
  % to transitions from base cases at t-1 to base cases at t. This allows
  % us to write a recursion. If the general form for SI(t),SF(t) is
  %     SI(t) = Li(t)*S0 + Mg(t)*x + Mh(t)*x,
  %     SF(t) = Lf(t)*S0 + Ng(t)*x + Nh(t)*x,
  % it turns out that the recursion is
  %     L(1) = D(1); Mg(1) = Mh(1) = 0; Ng(1) = G(1); Nh(1) = H(1);
  %     for t=2:nt
  %         Li(t) = D(t)*Lf(t-1) = D(t)*B1(t-1)*Li(t-1);
  %         Lf(t) = B1(t)*Li(t) = B1(t)*D(t)*Lf(t-1);
  %         Mg(t) = D(t)*Ng(t-1);
  %         Mh(t) = D(t)*Nh(t-1);
  %         Ng(t) = B1(t)*Mg(t) + B2(t)*G(t);
  %         Nh(t) = B1(t)*Mh(t) + B2(t)*H(t);
  %     end
  %
  % If SI,SF are organized first by blocks for each storage unit and within
  % the blocks by scenario, then the D matrix is simply made up by
  % repeating the str.tstep(t).TransMat matrix ns times in the diagonal
  % and then the columns weighted by the probabilities of the basecases at
  % t-1 given that we remained in basecases, and the rows are weighted by
  % the inverse of the probabilites of the scenarios at the beginning of t.
  % D(1) is special though, where each block is an nj(1) x 1 vector of ones.
  % The B matrices are formed by stacking appropriately sized diagonal
  % matrices for each storage unit along the diagonal, where each component
  % is simply the i-th element of beta times an identity matrix.
  if verbose
    fprintf('- Building expected storage-tracking mechanism.\n');
  end
  if ns
    % The following code assumes that no more variables will be added
    vv = om.get_idx();
    for t = 1:nt
      nsxnjt = ns*mdi.idx.nj(t);
      % Form G(t), H(t), B1(t), B2(t)
      G = sparse(nsxnjt, nvars);
      H = sparse(nsxnjt, nvars);
      B1 = sparse(nsxnjt, nsxnjt);
      B2 = sparse(nsxnjt, nsxnjt);
      for j = 1:mdi.idx.nj(t)
        ii  = ((1:ns)'-1)*mdi.idx.nj(t)+j;
        jj1 = (vv.i1.Psc(t,j,1):vv.iN.Psc(t,j,1))';
        jj2 = (vv.i1.Psd(t,j,1):vv.iN.Psd(t,j,1))';
        G = G + sparse(ii, jj1, -mdi.Delta_T  *  InEff(:,t), nsxnjt, nvars);
        H = H + sparse(ii, jj2, -mdi.Delta_T ./ OutEff(:,t), nsxnjt, nvars);
        B1 = B1 + sparse(ii, ii, beta1(:,t), nsxnjt, nsxnjt);
        B2 = B2 + sparse(ii, ii, beta2(:,t), nsxnjt, nsxnjt);
      end
      if t == 1
        % form Li, Lf, Mg, Mh, Ng, Nh, B1, B2 for t == 1
        jlist = [];
        for i=1:ns
          jlist = [ jlist; i*ones(mdi.idx.nj(t),1) ];
        end
        mdi.tstep(t).Li  = sparse((1:nsxnjt)', jlist, 1, nsxnjt, ns);
        mdi.tstep(t).Lf  = B1 * mdi.tstep(t).Li;
        mdi.tstep(t).Mg  = sparse(nsxnjt, nvars);   % Initial one is all zeros
        mdi.tstep(t).Mh  = sparse(nsxnjt, nvars);   % Initial one is all zeros
        mdi.tstep(t).Ng  = B2 * G;
        mdi.tstep(t).Nh  = B2 * H;
      else
        % Form D(t)
        D = sparse(nsxnjt, ns*mdi.idx.nj(t-1));
        p1 = mdi.CostWeights(1,1:mdi.idx.nj(t-1),t-1)';
        p1 = p1 / sum(p1);      % sigma(t)
        p2 = mdi.tstep(t).TransMat * p1;
        Di = spdiags(1./p2, 0, mdi.idx.nj(t), mdi.idx.nj(t)) * ...
                sparse(mdi.tstep(t).TransMat) * ...
                spdiags(p1, 0, mdi.idx.nj(t-1), mdi.idx.nj(t-1));
        for i = 1:ns
          D((i-1)*mdi.idx.nj(t)+1:i*mdi.idx.nj(t), (i-1)*mdi.idx.nj(t-1)+1:i*mdi.idx.nj(t-1)) = Di;
        end
        % Apply recursion, form Li, Lf, Mg, Mh, Ng, Nh
        mdi.tstep(t).Li = D  * mdi.tstep(t-1).Lf;
        mdi.tstep(t).Lf = B1 * mdi.tstep(t).Li;
        mdi.tstep(t).Mg = D * mdi.tstep(t-1).Ng;
        mdi.tstep(t).Mh = D * mdi.tstep(t-1).Nh;
        mdi.tstep(t).Ng = B1 * mdi.tstep(t).Mg + B2 * G;
        mdi.tstep(t).Nh = B1 * mdi.tstep(t).Mh + B2 * H;
      end
      mdi.tstep(t).G = G;
      mdi.tstep(t).H = H;
    end
  end

  % Now for the constraint indexing and creation.
  if verbose
    fprintf('- Building constraint submatrices.\n');
  end
  baseMVA = mdi.mpc.baseMVA;
  om.init_indexed_name('lin', 'Pmis', {nt, nj_max, nc_max+1});
  if mdi.DCMODEL
    % Construct all load flow equations using a DC flow model
    if verbose
      fprintf('  - Building DC flow constraints.\n');
    end
    om.init_indexed_name('lin', 'Pf', {nt, nj_max, nc_max+1});
    for t = 1:nt
      for j = 1:mdi.idx.nj(t)
        for k = 1:mdi.idx.nc(t,j)+1
          % First the flow constraints
          mpc = mdi.flow(t,j,k).mpc;
          ion = find(mpc.branch(:, BR_STATUS));
          [Bdc, Bl, Psh, PLsh] = makeBdc(baseMVA, mpc.bus, mpc.branch(ion,:));
          mdi.flow(t,j,k).PLsh = PLsh;     %% save for computing flows later
          negCg = sparse(mpc.gen(:,GEN_BUS), (1:ng)', -1, ...
                        mdi.idx.nb(t,j,k), ng);
          A = [Bdc negCg];
          b = -(mpc.bus(:,PD)+mpc.bus(:,GS))/baseMVA-Psh;
          vs = struct('name', {'Va', 'Pg'}, 'idx', {{t,j,k}, {t,j,k}});
          om.add_lin_constraint('Pmis', {t,j,k}, A, b, b, vs);
          % Then the thermal limits
          tmp = mpc.branch(ion,RATE_A)/baseMVA;
          iuncon = find(~tmp);
          tmp(iuncon) = Inf(size(iuncon));
          vs = struct('name', {'Va'}, 'idx', {{t,j,k}});
          om.add_lin_constraint('Pf', {t,j,k}, Bl, -tmp-PLsh, tmp-PLsh, vs);
        end
      end
    end
  else
    if verbose
      fprintf('  - Building load balance constraints.\n');
    end
    % Set simple generation - demand = 0 equations, one for each flow
    for t = 1:nt
      for j = 1:mdi.idx.nj(t)
        for k = 1:mdi.idx.nc(t,j)+1
          mpc = mdi.flow(t,j,k).mpc;
          A = sparse(ones(1, ng));
          b = 1.0*sum(mpc.bus(:, PD)+mpc.bus(:,GS))/baseMVA;
          vs = struct('name', {'Pg'}, 'idx', {{t,j,k}});
          om.add_lin_constraint('Pmis', {t,j,k}, A, b, b, vs);
        end
      end
    end
  end
  if mdi.IncludeFixedReserves
    if verbose
      fprintf('  - Building fixed zonal reserve constraints.\n');
    end
    om.init_indexed_name('lin', 'Pg_plus_R', {nt, nj_max, nc_max+1});
    om.init_indexed_name('lin', 'Rreq', {nt, nj_max, nc_max+1});
    for t = 1:nt
      for j = 1:mdi.idx.nj(t)
        for k = 1:mdi.idx.nc(t,j)+1
          % First the flow constraints
          mpc = mdi.flow(t,j,k).mpc;
          r = mdi.FixedReserves(t,j,k);
          ngr = length(r.igr);
          I = speye(ngr);
          Ar = sparse(1:ngr, r.igr, 1, ngr, ng);
          if UC
            A = [Ar I  ...
                  sparse(1:ngr, r.igr, -mpc.gen(r.igr, PMAX) / baseMVA, ngr, ng)];
            u = zeros(ngr, 1);
            vs = struct('name', {'Pg', 'R', 'u'}, 'idx', {{t,j,k}, {t,j,k}, {t}});
          else
            A = [Ar I];
            u = mpc.gen(r.igr, PMAX) / baseMVA;
            vs = struct('name', {'Pg', 'R'}, 'idx', {{t,j,k}, {t,j,k}});
          end
          om.add_lin_constraint('Pg_plus_R', {t,j,k}, A, [], u, vs);
          A = r.zones(:, r.igr);
          l = r.req / mpc.baseMVA;
          vs = struct('name', {'R'}, 'idx', {{t,j,k}});
          om.add_lin_constraint('Rreq', {t,j,k}, A, l, [], vs);
        end
      end
    end
  end
  
  % Set relationships between generator injections and charge/discharge
  % variables (-pg + psc + psd = 0)
  if verbose && ~isempty(mdi.Storage.UnitIdx)
    fprintf('  - Splitting storage injections into charge/discharge.\n');
  end
  om.init_indexed_name('lin', 'Ps', {nt, nj_max, nc_max+1});
  for t = 1:nt
    for j = 1:mdi.idx.nj(t)
      for k = 1:mdi.idx.nc(t,j)+1
        A = [sparse((1:ns)', mdi.Storage.UnitIdx, -1, ns, ng) Ins Ins];
        b = zeros(ns, 1);
        vs = struct('name', {'Pg', 'Psc', 'Psd'}, 'idx', {{t,j,k}, {t,j,k}, {t,j,k}});
        om.add_lin_constraint('Ps', {t,j,k}, A, b, b, vs);
      end
    end
  end
  
  % Construct y-variable restrictions on piecewise-linear costs. Note that
  % the restriction lines are computed using the full non-scaled cost in
  % gencost; any weighting of the cost must be then specified later in the
  % cost coefficients hitting the y variables (not 1 anymore).  Do it for
  % a complete c3sopf cell taking advantage of the fact that all p
  % injections are contiguous, as are all y variables for a c3sopf cell.
  % Also, note that makeAy assumes that every gen is online, which is
  % consistent with our formulation for unit commitment
  if verbose
    fprintf('  - Building CCV constraints for piecewise-linear costs.\n');
  end
  om.init_indexed_name('lin', 'ycon', {nt, nj_max, nc_max+1});
  for t = 1:nt,
    for j = 1:mdi.idx.nj(t)
      for k = 1:mdi.idx.nc(t,j)+1
        mpc = mdi.flow(t,j,k).mpc;
        [A, u] = makeAy(baseMVA, ng, mpc.gencost, 1, [], ng+1);
        vs = struct('name', {'Pg', 'y'}, 'idx', {{t,j,k}, {t,j,k}});
        om.add_lin_constraint('ycon', {t,j,k}, A, [], u, vs);
      end
    end
  end

  % The actual deviations from base flow must not exceed physical ramp rates
  % we'll get negative multiplier for right bound, fix when picking up
  % lambdas.
  % At issue: generators ousted in a contingency clearly violate this
  % transition; do not include constraints for these or for generators
  % whose commitment key is -1; either of these two possibilities will
  % result in a GEN_STATUS of 0, so we use that as the indicator.
  if verbose
    fprintf('  - Building contingency reserve constraints.\n');
  end
  om.init_indexed_name('lin', 'rampcont', {nt, nj_max, nc_max+1});
  for t =1:nt
    for j = 1:mdi.idx.nj(t)
      for k = 2:mdi.idx.nc(t,j)+1
        ii = find(mdi.flow(t,j,k).mpc.gen(:, GEN_STATUS) > 0);
        ngtmp = length(ii);
        A = sparse((1:ngtmp)', ii, 1, ngtmp, ng);
        u = mdi.flow(t,j,k).mpc.gen(ii,RAMP_10)/baseMVA;
        vs = struct('name', {'Pg', 'Pg'}, 'idx', {{t,j,1}, {t,j,k}});
        om.add_lin_constraint('rampcont', {t,j,k}, [-A A], -u, u, vs);
      end
    end
  end
  % The usual alpha-controlled equality of P0 and Pc does not make
  % sense when having many scenarios and hence many P0's .  Ditch.
  %
  % Positive reserve variables are larger than all increment variables in
  % all scenarios and flows of a given time slice 0 <= rpp - dpp; these
  % are the ones that set the price of reserves. Include all units that are
  % potentially committed.
  om.init_indexed_name('lin', 'dPpRp', {nt, nj_max, nc_max+1});
  for t = 1:nt
    for j = 1:mdi.idx.nj(t);
      for k = 1:mdi.idx.nc(t,j)+1
        ii = find(mdi.flow(t,j,k).mpc.gen(:, GEN_STATUS) > 0);
        ngtmp = length(ii);
        A = sparse((1:ngtmp)', ii, 1, ngtmp, ng);
        l = zeros(ngtmp, 1);
        vs = struct('name', {'dPp', 'Rpp'}, 'idx', {{t,j,k}, {t}});
        om.add_lin_constraint('dPpRp', {t,j,k}, [-A A], l, [], vs);
      end
    end
  end
  % Negative reserve variables are larger than all decrement variables in
  % all scenarios and flows of a given time slice  0 <= rpm - dpm; these
  % are the ones that set the price of reserves. Include all units that are
  % potentially committed.
  om.init_indexed_name('lin', 'dPmRm', {nt, nj_max, nc_max+1});
  for t = 1:nt
    for j = 1:mdi.idx.nj(t);
      for k = 1:mdi.idx.nc(t,j)+1
        ii = find(mdi.flow(t,j,k).mpc.gen(:, GEN_STATUS) > 0);
        ngtmp = length(ii);
        A = sparse((1:ngtmp)', ii, 1, ngtmp, ng);
        l = zeros(ngtmp, 1);
        vs = struct('name', {'dPm', 'Rpm'}, 'idx', {{t,j,k}, {t}});
        om.add_lin_constraint('dPmRm', {t,j,k}, [-A A], l, [], vs);
      end
    end
  end
  % The difference between the injection and the contract
  % is equal to the inc minus the dec: Ptjk - Ptc = dPp - dPm
  % Include all units that are potentially committed.
  om.init_indexed_name('lin', 'dPdef', {nt, nj_max, nc_max+1});
  for t = 1:nt
    for j = 1:mdi.idx.nj(t);
      for k = 1:mdi.idx.nc(t,j)+1
        ii = find(mdi.flow(t,j,k).mpc.gen(:, GEN_STATUS) > 0);
        ngtmp = length(ii);
        A = sparse((1:ngtmp)', ii, 1, ngtmp, ng);
        b = zeros(ngtmp, 1);
        vs = struct('name', {'Pg', 'Pc', 'dPp', 'dPm'}, ...
                    'idx', {{t,j,k}, {t}, {t,j,k}, {t,j,k}});
        om.add_lin_constraint('dPdef', {t,j,k}, [A -A -A A], b, b, vs);
      end
    end
  end
  
  % Go on to load following ramping restrictions.  Note that these
  % restrictions apply even if there is a change in the commitment status
  % of the generator.
  %
  % First, bound upward ramping reserves from below by all base-case
  % ramping possibilities, 0 <= rrp(t) -p(t+1)(j2)0 + p(t)(j1)0.  A ramping
  % reserve is needed at time t to be able to change to the needed dispatch at
  % time t+1.  An initial ramping reserve (t=0) would be data, not a
  % variable, and while ramping transitions from t=0 to t=1 should be
  % enforced, we do not allocate the reserve for t = 0.
  % Note: in the event that some future reserves are already locked in, we may
  %       not want to start at t = 1
  if verbose
    fprintf('  - Building ramping transitions and reserve constraints.\n');
  end
  om.init_indexed_name('lin', 'Rrp', {nt, nj_max, nj_max});
  % First, do from t=1:nt-1, since the last one is different and may not
  % even exist depending on the type of horizon
  for t = 1:nt-1
    for j1 = 1:mdi.idx.nj(t)      % j1 is at time t
      for j2 = 1:mdi.idx.nj(t+1)  % j2 is at time t+1
        if mdi.tstep(t+1).TransMask(j2,j1)
          A = [Ing -Ing Ing];
          l = zeros(ng, 1);
          vs = struct('name', {'Pg', 'Pg', 'Rrp'}, ...
                      'idx', {{t,j1,1}, {t+1,j2,1}, {t}});
          om.add_lin_constraint('Rrp', {t,j1,j2}, A, l, [], vs);
        end
      end
    end
  end
  % Now, pay special attention to a possible last type of ramping
  % constraint. If the horizon involves a terminal value at t=nt+1, then
  % this must also be enforced; in this case, additional ramping
  % reserves procured for t=nt must be defined.  If this
  % condition does not apply, then these reserves are not needed.
  if ~mdi.OpenEnded
    % pterminal <= rrp(nt) + p(nt,j1,0)
    for j1 = 1:mdi.idx.nj(nt)
      A = [Ing Ing];
      l = mdi.TerminalPg/baseMVA;
      vs = struct('name', {'Pg', 'Rrp'}, ...
                  'idx', {{nt,j1,1}, {nt}});
      om.add_lin_constraint('Rrp', {nt,j1,1}, A, l, [], vs);
    end
  end
  % Now on to downward ramping reserves.
  % First, bound downward ramping reserves from below by all base-case
  % ramping possibilities, 0 <= rrm(t) + p(t+1)j20 - p(t)j10
  om.init_indexed_name('lin', 'Rrm', {nt, nj_max, nj_max});
  % First, do from t=1:nt-1, since the last one is different and may not
  % even exist depending on the type of horizon
  for t = 1:nt-1
    for j1 = 1:mdi.idx.nj(t)      % j1 is at time t
      for j2 = 1:mdi.idx.nj(t+1)  % j2 is at time t+1
        if mdi.tstep(t+1).TransMask(j2,j1)
          A = [-Ing Ing Ing];
          l = zeros(ng, 1);
          vs = struct('name', {'Pg', 'Pg', 'Rrm'}, ...
                      'idx', {{t,j1,1}, {t+1,j2,1}, {t}});
          om.add_lin_constraint('Rrm', {t,j1,j2}, A, l, [], vs);
        end
      end
    end
  end
  % Now, pay special attention to a possible last type of ramping
  % constraint. If the horizon involves a terminal value at t=nt+1, then
  % this must also be enforced; in this case, additional ramping
  % reserves procured for t=nt must be defined.  If this
  % condition does not apply, then these reserves are not needed.
  if ~mdi.OpenEnded
    % -pterminal <= rrm(nt) - p(nt,j1,0)
    for j1 = 1:mdi.idx.nj(nt)
      A = [-Ing Ing];
      l = -mdi.TerminalPg/baseMVA;
      vs = struct('name', {'Pg', 'Rrm'}, ...
                  'idx', {{nt,j1,1}, {nt}});
      om.add_lin_constraint('Rrm', {nt,j1,1}, A, l, [], vs);
    end
  end
  %
  %
  % Now for the storage restrictions.
  if ns
    if verbose
      fprintf('  - Building storage constraints.\n');
    end
    % First bound sm(t) based on sm(t-1), with sm(1) being bound by the initial
    % data; this is for base case trajectories only
    om.init_indexed_name('lin', 'Sm', {nt, nj_max});
    if mdi.Storage.ForceCyclicStorage
      % sm(1) - beta1*s0 + beta2*Delta_T*[eta_c*psc(1,j,0) + (1/eta_d)*psd(1,j,0)] <= 0
      for j = 1:mdi.idx.nj(1)
        A = [ diagBeta2EtaIn1 diagBeta2overEtaOut1 Ins -spdiags(beta1(:,1), 0, ns, ns)];
        u = zeros(ns, 1);
        vs = struct('name', {'Psc', 'Psd', 'Sm', 'S0'}, 'idx', {{1,j,1}, {1,j,1}, {1}, {}});
        om.add_lin_constraint('Sm', {1,j}, A, [], u, vs);
      end
    else
      % sm(1) + beta2*Delta_T*[eta_c*psc(1,j,0) + (1/eta_d)*psd(1,j,0)] <= beta1*Initial/baseMVA
      for j = 1:mdi.idx.nj(1)
        A = [ diagBeta2EtaIn1 diagBeta2overEtaOut1 Ins ];
        u = beta1(:,1).*mdi.Storage.InitialStorageLowerBound/baseMVA;
        vs = struct('name', {'Psc', 'Psd', 'Sm'}, 'idx', {{1,j,1}, {1,j,1}, {1}});
        om.add_lin_constraint('Sm', {1,j}, A, [], u, vs);
      end
    end
    % Then the rest of the periods
    % sm(t) - beta1*(rho(t)*sm(t-1) + (1-rho(t))*s_I(t,j)) + beta2*Delta_T*[eta_c*psc(t,j,0) + (1/eta_d)*psd(t,j,0)] <= 0
    % where s_I(t,j) = L_I(t,j) * s0 + (Mg(t,j)+Mh(t,j)) * x
    for t = 2:nt
      for j = 1:mdi.idx.nj(t)
        Mj = mdi.tstep(t).Mg( j:mdi.idx.nj(t):(ns-1)*mdi.idx.nj(t)+j, :) + ...
             mdi.tstep(t).Mh( j:mdi.idx.nj(t):(ns-1)*mdi.idx.nj(t)+j, :);
        Lij = mdi.tstep(t).Li( j:mdi.idx.nj(t):(ns-1)*mdi.idx.nj(t)+j, :);
        diag1minusRhoBeta1 = spdiags((1-rho(:,t)) .* beta1(:,t), 0, ns, ns);
        A = sparse([1:ns,1:ns,1:ns,1:ns]', ...
                   [vv.i1.Psc(t,j,1):vv.iN.Psc(t,j,1), vv.i1.Psd(t,j,1):vv.iN.Psd(t,j,1), vv.i1.Sm(t-1):vv.iN.Sm(t-1), vv.i1.Sm(t):vv.iN.Sm(t)]', ...
                   [beta2EtaIn(:,t); beta2overEtaOut(:,t); -beta1(:,t).*rho(:,t); ones(ns,1)], ...
                   ns, nvars) ...
                - diag1minusRhoBeta1 * Mj;
        if mdi.Storage.ForceCyclicStorage
          As0 = sparse(ns, nvars);
          As0(:, vv.i1.S0:vv.iN.S0) = -diag1minusRhoBeta1 * Lij;
          A = A + As0;
          u = zeros(ns, 1);
        else
          u = full(diag1minusRhoBeta1 * Lij * mdi.Storage.InitialStorage/baseMVA);
        end
        om.add_lin_constraint('Sm', {t,j}, A, [], u);
      end
    end
    % Do the same we did for sm(t) for sp(t). First the initial step ...
    om.init_indexed_name('lin', 'Sp', {nt, nj_max});
    if mdi.Storage.ForceCyclicStorage
      % -sp(1) + beta1*s0 - beta2*Delta_T*[eta_c*psc(1,j,0) + (1/eta_d)*psd(1,j,0)] <= 0
      for j = 1:mdi.idx.nj(1)
        A = [ -diagBeta2EtaIn1 -diagBeta2overEtaOut1 -Ins spdiags(beta1(:,1), 0, ns, ns) ];
        u = zeros(ns, 1);
        vs = struct('name', {'Psc', 'Psd', 'Sp', 'S0'}, 'idx', {{1,j,1}, {1,j,1}, {1}, {}});
        om.add_lin_constraint('Sp', {1,j}, A, [], u, vs);
      end
    else
      % -sp(1) - beta2*Delta_T*[eta_c*psc(1,j,0) + (1/eta_d)*psd(1,j,0)] <= -beta1*Initial/baseMVA
      for j = 1:mdi.idx.nj(1)
        A = [ -diagBeta2EtaIn1 -diagBeta2overEtaOut1 -Ins ];
        u = -beta1(:,1).*mdi.Storage.InitialStorageUpperBound/baseMVA;
        vs = struct('name', {'Psc', 'Psd', 'Sp'}, 'idx', {{1,j,1}, {1,j,1}, {1}});
        om.add_lin_constraint('Sp', {1,j}, A, [], u, vs);
      end
    end
    % Then the rest of the periods
    % -sp(t) + beta1*(rho(t)*sp(t-1) + (1-rho(t))*s_I(t,j)) - beta2*Delta_T*[eta_c*psc(t,j,0) + (1/eta_d)*psd(t,j,0)] <= 0
    % where s_I(t,j) = L_I(t,j) * s0 + (Mg(t,j)+Mh(t,j)) * x
    for t = 2:nt
      for j = 1:mdi.idx.nj(t)
        Mj = mdi.tstep(t).Mg( j:mdi.idx.nj(t):(ns-1)*mdi.idx.nj(t)+j, :) + ...
             mdi.tstep(t).Mh( j:mdi.idx.nj(t):(ns-1)*mdi.idx.nj(t)+j, :);
        Lij = mdi.tstep(t).Li( j:mdi.idx.nj(t):(ns-1)*mdi.idx.nj(t)+j, :);
        diag1minusRhoBeta1 = spdiags((1-rho(:,t)) .* beta1(:,t), 0, ns, ns);
        A = sparse([1:ns,1:ns,1:ns,1:ns]', ...
                   [vv.i1.Psc(t,j,1):vv.iN.Psc(t,j,1), vv.i1.Psd(t,j,1):vv.iN.Psd(t,j,1), vv.i1.Sp(t-1):vv.iN.Sp(t-1), vv.i1.Sp(t):vv.iN.Sp(t)]', ...
                   [-beta2EtaIn(:,t); -beta2overEtaOut(:,t); beta1(:,t).*rho(:,t); -ones(ns,1)], ...
                   ns, nvars) ...
                + diag1minusRhoBeta1 * Mj;
        if mdi.Storage.ForceCyclicStorage
          As0 = sparse(ns, nvars);
          As0(:, vv.i1.S0:vv.iN.S0) = diag1minusRhoBeta1 * Lij;
          A = A + As0;
          u = zeros(ns, 1);
        else
          u = full(-diag1minusRhoBeta1 * Lij * mdi.Storage.InitialStorage/baseMVA);
        end
        om.add_lin_constraint('Sp', {t,j}, A, [], u);
      end
    end
    % Now go on and limit the amount of energy that can be used if a
    % contingency does happen. Bound sm first. First examine time period 1 wrt to initial
    % stored energy, then t=2 and on.
    om.init_indexed_name('lin', 'contSm', {nt, nj_max, nc_max+1});
    for j = 1:mdi.idx.nj(1)
      for k = 2:mdi.idx.nc(1,j)+1  %% NOTE NO k=1!!!
        if mdi.Storage.ForceCyclicStorage
          % beta4*Delta_T*[eta_c*psc(1,j,0) + (1/eta_d)*psd(1,j,0)] + beta3*Delta_T*[eta_c*psc(1,j,k) + (1/eta_d)*psd(1,j,k)] - beta5*s0 <= -sm_min(1)
          A = [ diagBeta4EtaIn1 diagBeta4overEtaOut1 diagBeta3EtaIn1 diagBeta3overEtaOut1 -spdiags(beta5(:,1), 0, ns, ns) ];
          u = -MinStorageLevel(:,1)/baseMVA;
          vs = struct('name', {'Psc', 'Psd', 'Psc', 'Psd', 'S0'}, ...
                      'idx', {{1,j,1}, {1,j,1}, {1,j,k}, {1,j,k}, {}});
        else
          % beta4*Delta_T*[eta_c*psc(1,j,0) + (1/eta_d)*psd(1,j,0)] + beta3*Delta_T*[eta_c*psc(1,j,k) + (1/eta_d)*psd(1,j,k)] <= beta5*Initial/baseMVA - sm_min(1)
          A = [ diagBeta4EtaIn1 diagBeta4overEtaOut1 diagBeta3EtaIn1 diagBeta3overEtaOut1 ];
          u = (beta5(:,1).*mdi.Storage.InitialStorageLowerBound - MinStorageLevel(:,1))/baseMVA;
          vs = struct('name', {'Psc', 'Psd', 'Psc', 'Psd'}, ...
                      'idx', {{1,j,1}, {1,j,1}, {1,j,k}, {1,j,k}});
        end
        om.add_lin_constraint('contSm', {1,j,k}, A, [], u, vs);
      end
    end
    % then the rest of the periods
    % -beta5*(rho(t)*sm(t-1) + (1-rho(t))*s_I(t,j)) + beta4*Delta_T*[eta_c*psc(t,j,0) + (1/eta_d)*psd(t,j,0)] + beta3*Delta_T*[eta_c*psc(t,j,k) + (1/eta_d)*psd(t,j,k)] <= -sm_min(t)
    % where s_I(t,j) = L_I(t,j) * s0 + (Mg(t,j)+Mh(t,j)) * x
    for t = 2:nt
      for j = 1:mdi.idx.nj(t)
        for k = 2:mdi.idx.nc(t,j)+1
          Mj = mdi.tstep(t).Mg( j:mdi.idx.nj(t):(ns-1)*mdi.idx.nj(t)+j, :) + ...
               mdi.tstep(t).Mh( j:mdi.idx.nj(t):(ns-1)*mdi.idx.nj(t)+j, :);
          Lij = mdi.tstep(t).Li( j:mdi.idx.nj(t):(ns-1)*mdi.idx.nj(t)+j, :);
          diag1minusRhoBeta5 = spdiags((1-rho(:,t)) .* beta5(:,t), 0, ns, ns);
          A = sparse([1:ns,1:ns,1:ns,1:ns,1:ns]', ...
                     [vv.i1.Psc(t,j,1):vv.iN.Psc(t,j,1), vv.i1.Psd(t,j,1):vv.iN.Psd(t,j,1), vv.i1.Psc(t,j,k):vv.iN.Psc(t,j,k), vv.i1.Psd(t,j,k):vv.iN.Psd(t,j,k), vv.i1.Sm(t-1):vv.iN.Sm(t-1)]', ...
                     [beta4EtaIn(:,t); beta4overEtaOut(:,t); beta3EtaIn(:,t); beta3overEtaOut(:,t); -beta5(:,t).*rho(:,t)], ...
                     ns, nvars) ...
                  - diag1minusRhoBeta5 * Mj;
          u = -MinStorageLevel(:,t)/baseMVA;
          if mdi.Storage.ForceCyclicStorage
            As0 = sparse(ns, nvars);
            As0(:, vv.i1.S0:vv.iN.S0) = -diag1minusRhoBeta5 * Lij;
            A = A + As0;
          else
            u = u + diag1minusRhoBeta5 * Lij * mdi.Storage.InitialStorageLowerBound/baseMVA;
          end
          om.add_lin_constraint('contSm', {t,j,k}, A, [], u);
        end
      end
    end
    % Bound sp first. First examine time period 1 wrt to initial
    % stored energy, then t=2 and on.
    om.init_indexed_name('lin', 'contSp', {nt, nj_max, nc_max+1});
    for j = 1:mdi.idx.nj(1)
      for k = 2:mdi.idx.nc(1,j)+1
        if mdi.Storage.ForceCyclicStorage
          % -beta4*Delta_T*[eta_c*psc(1,j,0) + (1/eta_d)*psd(1,j,0)] - beta3*Delta_T*[eta_c*psc(1,j,k) + (1/eta_d)*psd(1,j,k)] + beta5*s0 <= sp_max(1)
          A = [ -diagBeta4EtaIn1 -diagBeta4overEtaOut1 -diagBeta3EtaIn1 -diagBeta3overEtaOut1 spdiags(beta5(:,1), 0, ns, ns)];
          u = MaxStorageLevel(:,1)/baseMVA;
          vs = struct('name', {'Psc', 'Psd', 'Psc', 'Psd', 'S0'}, ...
                      'idx', {{1,j,1}, {1,j,1}, {1,j,k}, {1,j,k}, {}});
        else
          % -beta4*Delta_T*[eta_c*psc(1,j,0) + (1/eta_d)*psd(1,j,0)] - beta3*Delta_T*[eta_c*psc(1,j,k) + (1/eta_d)*psd(1,j,k)] <= -beta5*Initial/baseMVA + sp_max(1)
          A = [ -diagBeta4EtaIn1 -diagBeta4overEtaOut1 -diagBeta3EtaIn1 -diagBeta3overEtaOut1 ];
          u = (MaxStorageLevel(:,1) - beta5(:,1).*mdi.Storage.InitialStorageUpperBound)/baseMVA;
          vs = struct('name', {'Psc', 'Psd', 'Psc', 'Psd'}, ...
                      'idx', {{1,j,1}, {1,j,1}, {1,j,k}, {1,j,k}});
        end
        om.add_lin_constraint('contSp', {1,j,k}, A, [], u, vs);
      end
    end
    % then the rest of the periods
    % beta5*(rho(t)*sp(t-1) + (1-rho(t))*s_I(t,j)) - beta4*Delta_T*[eta_c*psc(t,j,0) + (1/eta_d)*psd(t,j,0)] - beta3*Delta_T*[eta_c*psc(t,j,k) + (1/eta_d)*psd(t,j,k)] <= sp_max(t)
    % where s_I(t,j) = L_I(t,j) * s0 + (Mg(t,j)+Mh(t,j)) * x
    for t = 2:nt
      for j = 1:mdi.idx.nj(t)
        for k = 2:mdi.idx.nc(t,j)+1
          Mj = mdi.tstep(t).Mg( j:mdi.idx.nj(t):(ns-1)*mdi.idx.nj(t)+j, :) + ...
               mdi.tstep(t).Mh( j:mdi.idx.nj(t):(ns-1)*mdi.idx.nj(t)+j, :);
          Lij = mdi.tstep(t).Li( j:mdi.idx.nj(t):(ns-1)*mdi.idx.nj(t)+j, :);
          diag1minusRhoBeta5 = spdiags((1-rho(:,t)) .* beta5(:,t), 0, ns, ns);
          A = sparse([1:ns,1:ns,1:ns,1:ns,1:ns]', ...
                     [vv.i1.Psc(t,j,1):vv.iN.Psc(t,j,1), vv.i1.Psd(t,j,1):vv.iN.Psd(t,j,1), vv.i1.Psc(t,j,k):vv.iN.Psc(t,j,k), vv.i1.Psd(t,j,k):vv.iN.Psd(t,j,k), vv.i1.Sp(t-1):vv.iN.Sp(t-1)]', ...
                     [-beta4EtaIn(:,t); -beta4overEtaOut(:,t); -beta3EtaIn(:,t); -beta3overEtaOut(:,t); beta5(:,t).*rho(:,t)], ...
                     ns, nvars) ...
                  + diag1minusRhoBeta5 * Mj;
          u = MaxStorageLevel(:,t)/baseMVA;
          if mdi.Storage.ForceCyclicStorage
            As0 = sparse(ns, nvars);
            As0(:, vv.i1.S0:vv.iN.S0) = diag1minusRhoBeta5 * Lij;
            A = A + As0;
          else
            u = u - diag1minusRhoBeta5 * Lij * mdi.Storage.InitialStorageUpperBound/baseMVA;
          end
          om.add_lin_constraint('contSp', {t,j,k}, A, [], u);
        end
      end
    end
  end
  
  % Now, if required, constrain the expected terminal storage quantity; two
  % different ways:
  if mdi.Storage.ForceExpectedTerminalStorage && mdi.Storage.ForceCyclicStorage
    error('most: ForceExpectedTerminalStorage and ForceCyclicStorage cannot be simultaneously true.');
  end
  if ns
    % The following code assumes that no more variables will be added
    if mdi.Storage.ForceExpectedTerminalStorage
      % 1) Constrain the expected terminal storage to be some target value
      A = sparse(ns, nvars);
      b = zeros(ns, 1);
      for j = 1:mdi.idx.nj(nt)
        Ngj = mdi.tstep(nt).Ng( j:mdi.idx.nj(nt):(ns-1)*mdi.idx.nj(nt)+j, :);
        Nhj = mdi.tstep(nt).Nh( j:mdi.idx.nj(nt):(ns-1)*mdi.idx.nj(nt)+j, :);
        Lfj = mdi.tstep(nt).Lf( j:mdi.idx.nj(nt):(ns-1)*mdi.idx.nj(nt)+j, :);
        A = A + mdi.CostWeights(1,j,nt) * (Ngj + Nhj);
        b = b + mdi.CostWeights(1,j,nt) * (Lfj * mdi.Storage.InitialStorage) / baseMVA;
      end
      endprob = sum(mdi.CostWeights(1,1:mdi.idx.nj(nt),nt)');
      A = (1/endprob) * A;
      b = (1/endprob) * b;
      l = mdi.Storage.ExpectedTerminalStorageMin / baseMVA - b;
      u = mdi.Storage.ExpectedTerminalStorageMax / baseMVA - b;
      om.add_lin_constraint('ESnt', A, l, u);
    elseif mdi.Storage.ForceCyclicStorage
      % 2) Constrain the initial storage (a variable) to be the same as the final expected storage
      A = sparse(ns, nvars);
      for j = 1:mdi.idx.nj(nt)
        Ngj = mdi.tstep(nt).Ng( j:mdi.idx.nj(nt):(ns-1)*mdi.idx.nj(nt)+j, :);
        Nhj = mdi.tstep(nt).Nh( j:mdi.idx.nj(nt):(ns-1)*mdi.idx.nj(nt)+j, :);
        A = A + mdi.CostWeights(1,j,nt) * (Ngj + Nhj);
      end
      endprob = sum(mdi.CostWeights(1,1:mdi.idx.nj(nt),nt)');
      A = (1/endprob) * A;
      for j = 1:mdi.idx.nj(nt)
        Lfj = mdi.tstep(nt).Lf(j:mdi.idx.nj(nt):(ns-1)*mdi.idx.nj(nt)+j, :);
        A(:, vv.i1.S0:vv.iN.S0) = A(:, vv.i1.S0:vv.iN.S0) ...
                  + mdi.CostWeights(1,j,nt) * Lfj;
      end
      A(:, vv.i1.S0:vv.iN.S0) = (1/endprob) * A(:, vv.i1.S0:vv.iN.S0) - speye(ns);
      b = zeros(ns, 1);
      om.add_lin_constraint('ESnt', A, b, b);
    end
  end

  % Dynamical system contraints
  if nzds || nyds
    if verbose
      fprintf('  - Building dynamical system constraints.\n');
    end
    % Compute matrices that give the expected dispatch in time period t
    % given that we make it to that period, for all generators at once,
    % when multiplied by x, i.e. E[p(t)] = E(t) * x
    for t = 1:nt
      E = sparse(ng, nvars);
      for j = 1:mdi.idx.nj(t)
        for k = 1:mdi.idx.nc(t,j)+1;
          E = E + (mdi.CostWeightsAdj(k,j,t)/mdi.StepProb(t)) * sparse((1:ng)', (vv.i1.Pg(t,j,k):vv.iN.Pg(t,j,k))', 1, ng, nvars);
        end
      end
      mdi.tstep(t).E = E;
    end
  end

  % Form the dynamical system state equations and bound constraints on the
  % state vector
  if nzds
    om.init_indexed_name('lin', 'DSz', {ntds-1});
    b = zeros(nzds, 1);
    for t = 1:ntds-1
      if t <= nt  % We have p(t) available to drive the dynamical system up to t=nt
        % Form the constraint matrix, so B*E*x + A*z(t) - I*z(t+1) = 0
        A = mdi.dstep(t).B * mdi.tstep(t).E;
      else
        % The dynamical system horizon is longer than the injection planning
        % horizon and we don't know what p(t) is, but continue to drive the
        % dynamical system as if p(t) = 0 and perhaps take that into account
        % when setting Ymax, Ymin in this time window.  That is, A*z(t) - I*z(t+1) = 0
        A = sparse(nzds, nvars);
      end
      A(:, vv.i1.Z(t):vv.iN.Z(t)) = mdi.dstep(t).A;
      A(:, vv.i1.Z(t+1):vv.iN.Z(t+1)) = -speye(nzds);
      om.add_lin_constraint('DSz', {t}, A, b, b);
    end
  end
  
  % Form the output equations and their restrictions
  if nyds
    om.init_indexed_name('lin', 'DSy', {ntds});
    for t = 1:ntds
      if t <= nt
        A = mdi.dstep(t).D * mdi.tstep(t).E;
      else
        A = sparse(nyds, nvars);
      end
      if nzds
        A(:, vv.i1.Z(t):vv.iN.Z(t)) = mdi.dstep(t).C;
      end
      l = mdi.dstep(t).ymin;
      u = mdi.dstep(t).ymax;
      om.add_lin_constraint('DSy', {t}, A, l, u);
    end
  end
  
  % UNIT COMMITMENT
  if UC
    if verbose
      fprintf('  - Building unit commitment constraints.\n');
    end
    % u(t,i) - u(t-1,i) - v(t,i) + w(t,i) = 0
    om.init_indexed_name('lin', 'uvw', {nt});
    for t = 1:nt
      if t == 1
        % First for t=1 when u(t-1,i) is really u(0,i) or u(nt,i)
        if mdi.UC.CyclicCommitment
          vs = struct('name', {'u', 'u', 'v', 'w'}, 'idx', {{1}, {nt}, {1}, {1}});
          A = [Ing -Ing -Ing Ing];
          b = zeros(ng, 1);
        else
          vs = struct('name', {'u', 'v', 'w'}, 'idx', {{1}, {1}, {1}});
          A = [Ing -Ing Ing];
          b = (mdi.UC.InitialState > 0);
        end
      else
        % Then for rest of periods
        vs = struct('name', {'u', 'u', 'v', 'w'}, 'idx', {{t}, {t-1}, {t}, {t}});
        A = [Ing -Ing -Ing Ing];
        b = zeros(ng, 1);
      end
      om.add_lin_constraint('uvw', {t}, A, b, b, vs);
    end
    % Then continue with minimimum up time constraints. Again, two
    % different forms depending on whether the horizon is cyclical or not
    om.init_indexed_name('lin', 'minup', {nt, ng});
    for t = 1:nt
      for i = 1:ng
        ti = t-mdi.UC.MinUp(i)+1:t;
        if mdi.UC.CyclicCommitment     % window is circular
          for tt = 1:length(ti)
            if ti(tt) < 1
              ti(tt) = nt + ti(tt);
            end
          end
        end
        % limit to positive time
        % even with CyclicCommitment, in case MinUp is longer than horizon
        % (which implies always ON or always OFF)
        ti = ti(ti>0);
        vs = struct('name', {'u'}, 'idx', {{t}});
        A = sparse(1, i, -1, 1, ng);
        for tt = 1:length(ti)
            vs(end+1).name = 'v';
            vs(end).idx  = {ti(tt)};
            A = [A sparse(1, i, 1, 1, ng)];
        end
        om.add_lin_constraint('minup', {t, i}, A, [], 0, vs);
      end
    end
    % Continue with minimimum downtime constraints. Two
    % different forms depending on whether the horizon is cyclical or not
    om.init_indexed_name('lin', 'mindown', {nt, ng});
    for t = 1:nt
      for i = 1:ng
        ti = t-mdi.UC.MinDown(i)+1:t;
        if mdi.UC.CyclicCommitment     % window is circular
          for tt = 1:length(ti)
            if ti(tt) < 1
              ti(tt) = nt + ti(tt);
            end
          end
        end
        % limit to positive time
        % even with CyclicCommitment, in case MinDown is longer than horizon
        % (which implies always ON or always OFF)
        ti = ti(ti>0);
        vs = struct('name', {'u'}, 'idx', {{t}});
        A = sparse(1, i, 1, 1, ng);
        for tt = 1:length(ti)
            vs(end+1).name = 'w';
            vs(end).idx  = {ti(tt)};
            A = [A sparse(1, i, 1, 1, ng)];
        end
        om.add_lin_constraint('mindown', {t, i}, A, [], 1, vs);
      end
    end
    % Limit generation ranges based on commitment status; first Pmax;
    % p - u*Pmax <= 0
    % For contingent flows, however, if a generator is ousted as a result
    % of the contingency, then this constraint should not be enforced.
    om.init_indexed_name('lin', 'uPmax', {nt, nj_max, nc_max+1});
    for t = 1:nt
      for j = 1:mdi.idx.nj(t)
        for k = 1:mdi.idx.nc(t,j)+1
          mpc = mdi.flow(t,j,k).mpc;
          ii = find(mpc.gen(:, GEN_STATUS));
          nii = length(ii);
          vs = struct('name', {'Pg', 'u'}, 'idx', {{t,j,k}, {t}});
          A = [ sparse(1:nii, ii, 1, nii, ng) ...
                sparse(1:nii, ii, -mpc.gen(ii, PMAX)/baseMVA, nii, ng) ];
          u = zeros(nii, 1);
          om.add_lin_constraint('uPmax', {t,j,k}, A, [], u, vs);
        end
      end
    end
    % Then Pmin,  -p + u*Pmin <= 0
    om.init_indexed_name('lin', 'uPmin', {nt, nj_max, nc_max+1});
    for t = 1:nt
      for j = 1:mdi.idx.nj(t)
        for k = 1:mdi.idx.nc(t,j)+1
          mpc = mdi.flow(t,j,k).mpc;
          ii = find(mpc.gen(:, GEN_STATUS));
          nii = length(ii);
          vs = struct('name', {'Pg', 'u'}, 'idx', {{t,j,k}, {t}});
          A = [ sparse(1:nii, ii, -1, nii, ng) ...
                sparse(1:nii, ii, mpc.gen(ii, PMIN)/baseMVA, nii, ng) ];
          u = zeros(nii, 1);
          om.add_lin_constraint('uPmin', {t,j,k}, A, [], u, vs);
        end
      end
    end
    % Then, if there is Qg coordination, do the same for Qg
    % q - u*Qmax <= 0
    % For contingent flows, however, if a generator is ousted as a result
    % of the contingency, then this constraint should not be enforced.
    if mdi.QCoordination
      om.init_indexed_name('lin', 'uQmax', {nt, nj_max, nc_max+1});
      for t = 1:nt
        for j = 1:mdi.idx.nj(t)
          for k = 1:mdi.idx.nc(t,j)+1
            mpc = mdi.flow(t,j,k).mpc;
            ii = find(mpc.gen(:, GEN_STATUS));
            nii = length(ii);
            vs = struct('name', {'Qg', 'u'}, 'idx', {{t,j,k}, {t}});
            A = [ sparse(1:nii, ii, 1, nii, ng) ...
                  sparse(1:nii, ii, -mpc.gen(ii, QMAX)/baseMVA, nii, ng) ];
            u = zeros(nii, 1);
            om.add_lin_constraint('uQmax', {t,j,k}, A, [], u, vs);
          end
        end
      end
      % Then Qmin,  -q + u*Qmin <= 0
      om.init_indexed_name('lin', 'uQmin', {nt, nj_max, nc_max+1});
      for t = 1:nt
        for j = 1:mdi.idx.nj(t)
          for k = 1:mdi.idx.nc(t,j)+1
            mpc = mdi.flow(t,j,k).mpc;
            ii = find(mpc.gen(:, GEN_STATUS));
            nii = length(ii);
            vs = struct('name', {'Qg', 'u'}, 'idx', {{t,j,k}, {t}});
            A = [ sparse(1:nii, ii, -1, nii, ng) ...
                  sparse(1:nii, ii, mpc.gen(ii, QMIN)/baseMVA, nii, ng) ];
            u = zeros(nii, 1);
            om.add_lin_constraint('uQmin', {t,j,k}, A, [], u, vs);
          end
        end
      end
    end
  end

  if verbose
    fprintf('- Building cost structures.\n');
  end
  % Start building the cost.  Two main components, the input data cost and
  % the coordination cost are specified.  The coordination cost is assumed to
  % have been buit with knowledge of the variable structure, and is simply
  % passed on.  The input data cost is assembled into the appropriate
  % spots.
  %
  % f = 0.5 * x' * (H1 + Hcoord) * x + (C1' + Ccoord) * x + c1 + ccoord

  % First assign the ramping costs; H1 has few coefficients initially and
  % this should make the shuffling and reordering of coefficients more
  % efficient.  All other accesses to H1 will be diagonal insertions, which
  % take less time than anti-diagonal insertions.
  % First do first period wrt to InitialPg.
  if mdi.OpenEnded
    om.init_indexed_name('qdc', 'RampWear', {nt, nj_max, nj_max});
  else
    om.init_indexed_name('qdc', 'RampWear', {nt+1, nj_max, nj_max});
  end
  for j = 1:mdi.idx.nj(1)
    w = mdi.tstep(1).TransMat(j,1);  % the probability of going from initial state to jth
    Q = spdiags(w * baseMVA^2 * mdi.RampWearCostCoeff(:,1), 0, ng, ng);
    c = -w * baseMVA * mdi.RampWearCostCoeff(:,1) .* mdi.InitialPg;
    vs = struct('name', {'Pg'}, 'idx', {{1,j,1}});
    k0 = w * 0.5 * mdi.RampWearCostCoeff(:,1)' * mdi.InitialPg.^2;
    om.add_quad_cost('RampWear', {1,j,1}, Q, c, k0, vs);
  end
  % Then the remaining periods
  for t = 2:nt
    for j2 = 1:mdi.idx.nj(t)
      for j1 = 1:mdi.idx.nj(t-1)
        w = mdi.tstep(t).TransMat(j2,j1) * mdi.CostWeights(1, j1, t-1);
        h = w * baseMVA^2 * mdi.RampWearCostCoeff(:,t);
        i = (1:ng)';
        j = ng+(1:ng)';
        Q = sparse([i;j;i;j], [i;i;j;j], [h;-h;-h;h], 2*ng, 2*ng);
        vs = struct('name', {'Pg', 'Pg'}, 'idx', {{t-1,j1,1}, {t,j2,1}});
        om.add_quad_cost('RampWear', {t,j1,j2}, Q, zeros(2*ng,1), 0, vs);
      end
    end
  end
  % Finally, if there is a terminal state problem, apply cost to
  % the transition starting from t=nt.  Note that in this case
  % mdi.tstep(nt+1).TransMat must be defined! it is the only piece of data
  % that makes sense for nt+1; all other fields in mdi.tstep(nt+1) can be empty.
  if ~mdi.OpenEnded
    for j = 1:mdi.idx.nj(nt)
      w = mdi.tstep(nt+1).TransMat(1, j) * mdi.CostWeights(1, j, nt);
      Q = spdiags(w * baseMVA^2 * mdi.RampWearCostCoeff(:,nt+1), 0, ng, ng);
      c = -w * baseMVA * mdi.RampWearCostCoeff(:,nt+1) .* mdi.TerminalPg;
      vs = struct('name', {'Pg'}, 'idx', {{nt,j,1}});
      k0 = w * 0.5 * mdi.RampWearCostCoeff(:,nt+1)' * mdi.TerminalPg.^2;
      om.add_quad_cost('RampWear', {nt+1,j,1}, Q, c, k0, vs);
    end
  end

  % Now go on and assign energy, inc/dec and contingency reserves
  % costs for all committed units.
  om.init_indexed_name('qdc', 'Cp', {nt, nj_max, nc_max+1});
  om.init_indexed_name('qdc', 'Cy', {nt, nj_max, nc_max+1});
  om.init_indexed_name('qdc', 'Cpp', {nt, nj_max, nc_max+1});
  om.init_indexed_name('qdc', 'Cpm', {nt, nj_max, nc_max+1});
  if mdi.IncludeFixedReserves
    om.init_indexed_name('qdc', 'Rcost', {nt, nj_max, nc_max+1});
  end
  om.init_indexed_name('qdc', 'Crpp', {nt});
  om.init_indexed_name('qdc', 'Crpm', {nt});
  for t = 1:nt
    for j = 1:mdi.idx.nj(t)
      for k = 1:mdi.idx.nc(t,j)+1
        w = mdi.CostWeightsAdj(k,j,t);     %% NOTE (k,j,t) order !!!

        % weighted polynomial energy costs for committed units
        gc = mdi.flow(t,j,k).mpc.gencost;
        ipol = find(gc(:, MODEL) == POLYNOMIAL);
        if ~isempty(ipol)
          ncost = gc(ipol(1), NCOST);
          if all(gc(ipol, NCOST) == ncost)    %% uniform order of polynomials
            %% use vectorized code
            if ncost > 3
              error('most: polynomial generator costs of order higher than quadratic not supported');
            elseif ncost == 3
              Q = sparse(ipol, ipol, 2 * w * baseMVA^2*gc(ipol, COST), ng, ng);
            else
              Q = sparse(ng,ng);
            end
            c = zeros(ng, 1);
            if ncost >= 2
              c(ipol) = w * baseMVA*gc(ipol, COST+ncost-2);
            end
            k0 = w * sum(gc(ipol, COST+ncost-1));
          else                                %% non-uniform order of polynomials
            %% use a loop
            Q = sparse(ng,ng);
            c = zeros(ng, 1);
            for i = ipol'
              ncost = gc(i, NCOST);
              if ncost > 3
                error('most: polynomial generator costs of order higher than quadratic not supported');
              elseif ncost == 3
                Q(i,i) = 2 * w * baseMVA^2*gc(i, COST);
              end
              if ncost >= 2
                c(i) = w * baseMVA*gc(i, COST+ncost-2);
              end
              k0 = w * gc(i, COST+ncost-1);
            end
          end
          vs = struct('name', {'Pg'}, 'idx', {{t,j,k}});
          om.add_quad_cost('Cp', {t,j,k}, Q, c, k0, vs);
        end

        % weighted y-variables for piecewise linear energy costs for committed units
        % ipwl = find( (mdi.flow(t,j,k).mpc.gen(:,GEN_STATUS) > 0) & (gc(:,MODEL) == PW_LINEAR));
        if mdi.idx.ny(t,j,k)
          c = w * ones(mdi.idx.ny(t,j,k),1);
          vs = struct('name', {'y'}, 'idx', {{t,j,k}});
          om.add_quad_cost('Cy', {t,j,k}, [], c, 0, vs);
        end

        % inc and dec offers for each flow
        c = w * baseMVA * mdi.offer(t).PositiveActiveDeltaPrice(:);
        vs = struct('name', {'dPp'}, 'idx', {{t,j,k}});
        om.add_quad_cost('Cpp', {t,j,k}, [], c, 0, vs);
        c = w * baseMVA * mdi.offer(t).NegativeActiveDeltaPrice(:);
        vs = struct('name', {'dPm'}, 'idx', {{t,j,k}});
        om.add_quad_cost('Cpm', {t,j,k}, [], c, 0, vs);

        % weighted fixed reserves cost
        if mdi.IncludeFixedReserves
          c = w * mdi.FixedReserves(t,j,k).cost(r.igr) * baseMVA;
          vs = struct('name', {'R'}, 'idx', {{t,j,k}});
          om.add_quad_cost('Rcost', {t,j,k}, [], c, 0, vs);
        end
      end
    end
    
    % contingency reserve costs
    c = baseMVA * mdi.StepProb(t) * mdi.offer(t).PositiveActiveReservePrice(:);
    vs = struct('name', {'Rpp'}, 'idx', {{t}});
    om.add_quad_cost('Crpp', {t}, [], c, 0, vs);
    c = baseMVA * mdi.StepProb(t) * mdi.offer(t).NegativeActiveReservePrice(:);
    vs = struct('name', {'Rpm'}, 'idx', {{t}});
    om.add_quad_cost('Crpm', {t}, [], c, 0, vs);
  end
  % Assign load following ramp reserve costs.  Do first nt-1 periods first
  om.init_indexed_name('qdc', 'Crrp', {mdi.idx.ntramp});
  om.init_indexed_name('qdc', 'Crrm', {mdi.idx.ntramp});
  for t = 1:nt-1,
    c = baseMVA * mdi.StepProb(t+1) * mdi.offer(t).PositiveLoadFollowReservePrice(:);
    vs = struct('name', {'Rrp'}, 'idx', {{t}});
    om.add_quad_cost('Crrp', {t}, [], c, 0, vs);
    c = baseMVA * mdi.StepProb(t+1) * mdi.offer(t).NegativeLoadFollowReservePrice(:);
    vs = struct('name', {'Rrm'}, 'idx', {{t}});
    om.add_quad_cost('Crrm', {t}, [], c, 0, vs);
  end
  % Then do last period if needed Terminal state case
  if ~mdi.OpenEnded
    %% are these costs missing a mdi.StepProb(t)?  -- rdz
    c = baseMVA * mdi.offer(nt).PositiveLoadFollowReservePrice(:);
    vs = struct('name', {'Rrp'}, 'idx', {{nt}});
    om.add_quad_cost('Crrp', {nt}, [], c, 0, vs);
    c = baseMVA * mdi.offer(nt).NegativeLoadFollowReservePrice(:);
    vs = struct('name', {'Rrm'}, 'idx', {{nt}});
    om.add_quad_cost('Crrm', {nt}, [], c, 0, vs);
  end
  % Assign startup/shutdown costs, if any, and fixed operating costs
  if UC
    om.init_indexed_name('qdc', 'c00', {nt});
    om.init_indexed_name('qdc', 'startup', {nt});
    om.init_indexed_name('qdc', 'shutdown', {nt});
    for t = 1:nt
      ww = zeros(ng, 1);
      for j = 1:mdi.idx.nj(t)
        for k = 1:mdi.idx.nc(t)+1
          ww = ww + mdi.CostWeightsAdj(k,j,t) * mdi.flow(t,j,k).mpc.gen(:, GEN_STATUS);
        end
      end
      c = ww.*mdi.UC.c00(:,t);
      vs = struct('name', {'u'}, 'idx', {{t}});
      om.add_quad_cost('c00', {t}, [], c, 0, vs);
      c = mdi.StepProb(t)*mdi.flow(t,1,1).mpc.gencost(:, STARTUP);
      vs = struct('name', {'v'}, 'idx', {{t}});
      om.add_quad_cost('startup', {t}, [], c, 0, vs);
      c = mdi.StepProb(t)*mdi.flow(t,1,1).mpc.gencost(:, SHUTDOWN);
      vs = struct('name', {'w'}, 'idx', {{t}});
      om.add_quad_cost('shutdown', {t}, [], c, 0, vs);
    end
  end
  % Finally, assign any value to leftover stored energy
  if ns
    A1 = sparse(ns, ns);
    A2 = sparse(ns, nvars);
    A3 = sparse(ns, nvars);
    A4 = sparse(ns, nvars);
    A5 = sparse(ns, nvars);
    A6 = sparse(ns, nvars);
    A7 = sparse(ns, nvars);
    % The following code assumes that no more variables will be added
    vv = om.get_idx();
    for t = 1:nt
      % Compute cost coefficients for value of expected leftover storage
      % after a contingency
      for j = 1:mdi.idx.nj(t)
        % pick rows for jth base injections
        Gtj0  = mdi.tstep(t).G(  j:mdi.idx.nj(t):(ns-1)*mdi.idx.nj(t)+j, :);
        Htj0  = mdi.tstep(t).H(  j:mdi.idx.nj(t):(ns-1)*mdi.idx.nj(t)+j, :);
        Litj0 = mdi.tstep(t).Li( j:mdi.idx.nj(t):(ns-1)*mdi.idx.nj(t)+j, :);
        Mgtj0 = mdi.tstep(t).Mg( j:mdi.idx.nj(t):(ns-1)*mdi.idx.nj(t)+j, :);
        Mhtj0 = mdi.tstep(t).Mh( j:mdi.idx.nj(t):(ns-1)*mdi.idx.nj(t)+j, :);
        sum_psi_tjk = sum(mdi.CostWeights(2:mdi.idx.nc(t,j)+1,j,t));
        if t == nt
          A1 = A1 + mdi.CostWeights(1,j,t) * spdiags(OutEff(:,t) .* beta1(:,t), 0, ns, ns) * Litj0;
          A2 = A2 + mdi.CostWeights(1,j,t) * spdiags(OutEff(:,t) .* beta1(:,t), 0, ns, ns) * Mgtj0;
          A3 = A3 + mdi.CostWeights(1,j,t) * spdiags(OutEff(:,t) .* beta1(:,t), 0, ns, ns) * Mhtj0;
          A4 = A4 + mdi.CostWeights(1,j,t) * spdiags(OutEff(:,t) .* beta2(:,t), 0, ns, ns) * Gtj0;
          A5 = A5 + mdi.CostWeights(1,j,t) * spdiags(OutEff(:,t) .* beta2(:,t), 0, ns, ns) * Htj0;
        end
        A1 = A1 + sum_psi_tjk * spdiags(OutEff(:,t) .* beta5(:,t), 0, ns, ns) * Litj0;
        A2 = A2 + sum_psi_tjk * (spdiags(OutEff(:,t) .* beta5(:,t), 0, ns, ns) * Mgtj0 + spdiags(OutEff(:,t) .* beta4(:,t), 0, ns, ns) * Gtj0);
        A3 = A3 + sum_psi_tjk * (spdiags(OutEff(:,t) .* beta5(:,t), 0, ns, ns) * Mhtj0 + spdiags(OutEff(:,t) .* beta4(:,t), 0, ns, ns) * Htj0);
        for k = 2:mdi.idx.nc(t,j)+1
          ii  = (1:ns)';
          jj1 = (vv.i1.Psc(t,j,k):vv.iN.Psc(t,j,k))';
          jj2 = (vv.i1.Psd(t,j,k):vv.iN.Psd(t,j,k))';
          Gtjk = sparse(ii, jj1, -mdi.Delta_T  *  InEff(:,t), ns, nvars);
          Htjk = sparse(ii, jj2, -mdi.Delta_T ./ OutEff(:,t), ns, nvars);
          A6 = A6 + mdi.CostWeights(k,j,t) * spdiags(OutEff(:,t) .* beta3(:,t), 0, ns, ns) * Gtjk;
          A7 = A7 + mdi.CostWeights(k,j,t) * spdiags(OutEff(:,t) .* beta3(:,t), 0, ns, ns) * Htjk;
        end
      end
    end
    Cfstor = -baseMVA * ...
       (mdi.Storage.TerminalStoragePrice'      * (A2 + A3) + ...
        mdi.Storage.TerminalChargingPrice0'    * A4 + ...
        mdi.Storage.TerminalDischargingPrice0' * A5 + ...
        mdi.Storage.TerminalChargingPriceK'    * A6 + ...
        mdi.Storage.TerminalDischargingPriceK' * A7);
    if mdi.Storage.ForceCyclicStorage
      % If the horizon model for the storage is cyclic and therefore s0 is a
      % variable, then that initial storage must come at a cost,
      % (InitialStorageCost) otherwise the optimizer will force the gratis
      % s0 up just to have (possibly) more storage left at the end.
      Cfstor(vv.i1.S0:vv.iN.S0) = ...
        Cfstor(vv.i1.S0:vv.iN.S0) + ...
            baseMVA * mdi.Storage.InitialStorageCost';
      % and the term in the final expected storage related to s0 is also
      % not constant, so must be included in the objective function
      Cfstor(vv.i1.S0:vv.iN.S0) = ...
          Cfstor(vv.i1.S0:vv.iN.S0) - ...
          baseMVA * mdi.Storage.TerminalStoragePrice' * A1;
    end
    om.add_quad_cost('fstor', [], Cfstor', 0);

    % The following is a hack to make the storage state bounds tight;
    % assign them a very small cost
    om.init_indexed_name('qdc', 'SpSmFudge', {nt});
    c = 1e-2 * [-ones(ns,1); ones(ns,1)];
    for t = 1:nt
      vs = struct('name', {'Sm', 'Sp'}, 'idx', {{t}, {t}});
      om.add_quad_cost('SpSmFudge', {t}, [], c, 0, vs);
    end
  else
    Cfstor = sparse(1, nvars);
  end

  % Plug into struct
  if verbose
    fprintf('- Assembling full set of costs.\n');
  end
  [Q, c, k0] = om.params_quad_cost();
  mdi.QP.Cfstor = Cfstor;
  mdi.QP.H1 = Q;
  mdi.QP.C1 = c;
  mdi.QP.c1 = k0;
end     % if mpopt.most.build_model

% With all pieces of the cost in place, can proceed to build the total
% cost now.
mdi.QP.H = mdi.QP.H1;
mdi.QP.C = mdi.QP.C1;
mdi.QP.c = mdi.QP.c1;
if isfield(mdi, 'CoordCost') && ...
        (~isempty(mdi.CoordCost.Cuser) || ~isempty(mdi.CoordCost.Huser))
  if verbose
    fprintf('- Adding coordination cost to standard cost.\n');
  end
  nvuser = length(mdi.CoordCost.Cuser);
  nvars = mdi.idx.nvars;
  mdi.QP.H = mdi.QP.H + ...
            [ mdi.CoordCost.Huser       sparse(nvuser,nvars-nvuser) ;
            sparse(nvars-nvuser,nvuser)  sparse(nvars-nvuser,nvars-nvuser) ];
  mdi.QP.C(1:nvuser) = mdi.QP.C(1:nvuser) +  mdi.CoordCost.Cuser(:);
  mdi.QP.c = mdi.QP.c + mdi.CoordCost.cuser;
  
%   cp = struct('Cw', mdi.CoordCost.Cuser(:), ...
%         'H', [ mdi.CoordCost.Huser     sparse(nvuser,nvars-nvuser) ;
%             sparse(nvars-nvuser,nvuser) sparse(nvars-nvuser,nvars-nvuser) ]);
%   om.add_legacy_cost('CoordCost', cp);
end

mdi.om = om;
[vv, ll] = om.get_idx();
if verbose
  fprintf('- Assembling full set of constraints.\n');
end
[mdi.QP.A, mdi.QP.l, mdi.QP.u] = om.params_lin_constraint();
if verbose
  fprintf('- Assembling full set of variable bounds.\n');
end
[mdi.QP.x0, mdi.QP.xmin, mdi.QP.xmax, mdi.QP.vtype] = om.params_var();

et_setup = toc(t0);
t0 = tic;

% Call solver!
mdo = mdi;
if mpopt.most.solve_model
  %% check consistency of model options (in case mdi was built in previous call)
  if mdi.DCMODEL ~= mo.DCMODEL
    error('MDI.DCMODEL inconsistent with MPOPT.most.dc_model');
  end
  if mdi.IncludeFixedReserves ~= mo.IncludeFixedReserves
    error('MDI.IncludeFixedReserves inconsistent with MPOPT.most.fixed_res (and possible presence of MDI.FixedReserves(t,j,k))');
  end
  if mdi.SecurityConstrained ~= mo.SecurityConstrained
    error('MDI.SecurityConstrained inconsistent with MPOPT.most.security_constraints (and possible presence of MDI.cont(t,j).contab)');
  end
  if mdi.QCoordination ~= mo.QCoordination
    error('MDI.QCoordination inconsistent with MPOPT.most.q_coordination');
  end
  if mdi.Storage.ForceCyclicStorage ~= mo.ForceCyclicStorage
    error('MDI.Storage.ForceCyclicStorage inconsistent with MPOPT.most.storage.cyclic');
  end
  if mdi.Storage.ForceExpectedTerminalStorage ~= mo.ForceExpectedTerminalStorage
    error('MDI.Storage.ForceExpectedTerminalStorage inconsistent with MPOPT.most.storage.terminal_target (and possible presence of MDI.Storage.ExpectedTerminalStorageAim|Min|Max)');
  end
  if mdi.UC.run ~= UC
    error('MDI.UC.run inconsistent with MPOPT.most.uc.run (and possible presence of MDI.UC.CommitKey)');
  end
  %% set options
  if any(any(mdi.QP.H))
    model = 'QP';
  else
    model = 'LP';
  end
  if UC
    model = ['MI' model];
  end
  mdo.QP.opt = mpopt2qpopt(mpopt, model, 'most');
  if verbose
    fprintf('- Calling %s solver.\n\n', model);
    fprintf('============================================================================\n\n');
  end
  if UC
    [mdo.QP.x, mdo.QP.f, mdo.QP.exitflag, mdo.QP.output, ...
            mdo.QP.lambda ] = miqps_matpower( mdi.QP.H, mdi.QP.C, ...
                mdi.QP.A, mdi.QP.l, mdi.QP.u, mdi.QP.xmin, mdi.QP.xmax, ...
                [], mdi.QP.vtype, mdo.QP.opt);
  else
    [mdo.QP.x, mdo.QP.f, mdo.QP.exitflag, mdo.QP.output, ...
            mdo.QP.lambda ] = qps_matpower( mdi.QP.H, mdi.QP.C, ...
                mdi.QP.A, mdi.QP.l, mdi.QP.u, mdi.QP.xmin, mdi.QP.xmax, ...
                [], mdo.QP.opt);
  end
  if mdo.QP.exitflag > 0
    success = 1;
    if verbose
      fprintf('\n============================================================================\n');
      fprintf('- MOST: %s solved successfully.\n', model);
    end
  else
    success = 0;
    if verbose
      fprintf('\n============================================================================\n');
      fprintf('- MOST: %s solver ''%s'' failed with exit flag = %d\n', model, mdo.QP.opt.alg, mdo.QP.exitflag);
    end
%     fprintf('\n============================================================================\n');
%     fprintf('- MOST: %s solver ''%s'' failed with exit flag = %d\n', model, mdo.QP.opt.alg, mdo.QP.exitflag);
%     fprintf('  You can query the workspace to debug.\n')
%     fprintf('  When finished, type the word "dbcont" to continue.\n\n');
%     keyboard;
  end
  % Unpack results
  if verbose
    fprintf('- Post-processing results.\n');
  end
  if success
    for t = 1:nt
      if UC
        mdo.UC.CommitSched(:, t) = mdo.QP.x(vv.i1.u(t):vv.iN.u(t));
      end
      for j = 1:mdi.idx.nj(t)
        for k = 1:mdi.idx.nc(t,j)+1
          mpc = mdo.flow(t,j,k).mpc;      %% pull mpc from output struct
          % Some initialization of data
          if mdo.DCMODEL
            mpc.bus(:, VM) = 1;
          end
          % Injections and shadow prices
          mpc.gen(:, PG) = baseMVA * mdo.QP.x(vv.i1.Pg(t,j,k):vv.iN.Pg(t,j,k));
          %% need to update Qg for loads consistent w/constant power factor
          Pmin = mpc.gen(:, PMIN);
          Qmin = mpc.gen(:, QMIN);
          Qmax = mpc.gen(:, QMAX);
          ivl = find( isload(mpc.gen) & (Qmin ~= 0 | Qmax ~= 0) );
          Qlim = (Qmin(ivl) == 0) .* Qmax(ivl) + (Qmax(ivl) == 0) .* Qmin(ivl);
          mpc.gen(ivl, QG) = mpc.gen(ivl, PG) .* Qlim ./ Pmin(ivl);
          if mdo.DCMODEL
            %% bus angles
            mpc.bus(:, VA) = (180/pi) * mdo.QP.x(vv.i1.Va(t,j,k):vv.iN.Va(t,j,k));
          
            %% nodal prices
            price = (mdo.QP.lambda.mu_u(ll.i1.Pmis(t,j,k):ll.iN.Pmis(t,j,k))-mdo.QP.lambda.mu_l(ll.i1.Pmis(t,j,k):ll.iN.Pmis(t,j,k))) / baseMVA;
            mpc.bus(:, LAM_P) = price;
          
            %% line flows and line limit shadow prices
            mpc.branch(:, PF) = 0;
            mpc.branch(:, QF) = 0;
            mpc.branch(:, PT) = 0;
            mpc.branch(:, QT) = 0;
            mpc.branch(:, MU_SF) = 0;
            mpc.branch(:, MU_ST) = 0;
            ion = find(mpc.branch(:, BR_STATUS));
            rows = ll.i1.Pf(t,j,k):ll.iN.Pf(t,j,k);
            cols = vv.i1.Va(t,j,k):vv.iN.Va(t,j,k);
            lf = baseMVA * (mdo.QP.A(rows,cols) * mdo.QP.x(cols) + mdo.flow(t,j,k).PLsh);
            mpc.branch(ion, PF) = lf;
            mpc.branch(ion, PT) = -lf;
            mpc.branch(ion, MU_SF) = mdo.QP.lambda.mu_u(rows) / baseMVA;
            mpc.branch(ion, MU_ST) = mdo.QP.lambda.mu_l(rows) / baseMVA;
          else
            %% system price
            price = (mdo.QP.lambda.mu_l(ll.i1.Pmis(t,j,k):ll.iN.Pmis(t,j,k))-mdo.QP.lambda.mu_u(ll.i1.Pmis(t,j,k):ll.iN.Pmis(t,j,k))) / baseMVA;
            mpc.bus(:, LAM_P) = price;
          end
          if UC
            % igenon does not contain gens ousted because of a contingency or
            % a forced-off UC.CommitKey
            igenon = find(mpc.gen(:, GEN_STATUS));
            u = mdo.QP.x(vv.i1.u(t):vv.iN.u(t));
            mpc.gen(igenon, GEN_STATUS) = u(igenon);
            gs = mpc.gen(igenon, GEN_STATUS) > 0; % gen status
            mpc.gen(:, MU_PMAX) = 0;
            mpc.gen(:, MU_PMIN) = 0;
            mpc.gen(igenon, MU_PMAX) = gs .* ...
                    mdo.QP.lambda.mu_u(ll.i1.uPmax(t,j,k):ll.iN.uPmax(t,j,k)) / baseMVA;
            mpc.gen(igenon, MU_PMIN) = gs .* ...
                    mdo.QP.lambda.mu_u(ll.i1.uPmin(t,j,k):ll.iN.uPmin(t,j,k)) / baseMVA;
            if mdo.QCoordination
              mpc.gen(:, MU_QMAX) = 0;
              mpc.gen(:, MU_QMIN) = 0;
              mpc.gen(igenon, MU_QMAX) = gs .* ...
                      mdo.QP.lambda.mu_u(ll.i1.uQmax(t,j,k):ll.iN.uQmax(t,j,k)) / baseMVA;
              mpc.gen(igenon, MU_QMIN) = gs .* ...
                      mdo.QP.lambda.mu_u(ll.i1.uQmin(t,j,k):ll.iN.uQmin(t,j,k)) / baseMVA;
            end
          else
            gs = mpc.gen(:, GEN_STATUS) > 0;      % gen status
            mpc.gen(:, MU_PMAX) = gs .* ...
                    mdo.QP.lambda.upper(vv.i1.Pg(t,j,k):vv.iN.Pg(t,j,k)) / baseMVA;
            mpc.gen(:, MU_PMIN) = gs .* ...
                    mdo.QP.lambda.lower(vv.i1.Pg(t,j,k):vv.iN.Pg(t,j,k)) / baseMVA;
            if mdo.QCoordination
              mpc.gen(:, MU_QMAX) = gs .* ...
                      mdo.QP.lambda.upper(vv.i1.Qg(t,j,k):vv.iN.Qg(t,j,k)) / baseMVA;
              mpc.gen(:, MU_QMIN) = gs .* ...
                      mdo.QP.lambda.lower(vv.i1.Qg(t,j,k):vv.iN.Qg(t,j,k)) / baseMVA;
            end
          end
          if mdi.IncludeFixedReserves
            z = zeros(ng, 1);
            r = mdo.FixedReserves(t,j,k);
            r.R   = z;
            r.prc = z;
            r.mu = struct('l', z, 'u', z, 'Pmax', z);
            r.totalcost = sum(om.eval_quad_cost(mdo.QP.x, 'Rcost', {t,j,k}));
            r.R(r.igr) = mdo.QP.x(vv.i1.R(t,j,k):vv.iN.R(t,j,k)) * baseMVA;
            for gg = r.igr
              iz = find(r.zones(:, gg));
              kk = ll.i1.Rreq(t,j,k):ll.iN.Rreq(t,j,k);
              r.prc(gg) = sum(mdo.QP.lambda.mu_l(kk(iz))) / baseMVA;
            end
            r.mu.l(r.igr)    = mdo.QP.lambda.lower(vv.i1.R(t,j,k):vv.iN.R(t,j,k)) / baseMVA;
            r.mu.u(r.igr)    = mdo.QP.lambda.upper(vv.i1.R(t,j,k):vv.iN.R(t,j,k)) / baseMVA;
            r.mu.Pmax(r.igr) = mdo.QP.lambda.mu_u(ll.i1.Pg_plus_R(t,j,k):ll.iN.Pg_plus_R(t,j,k)) / baseMVA;
            mpc.reserves = r;
          end
          mdo.flow(t,j,k).mpc = mpc;     %% stash modified mpc in output struct
        end
      end
      % Contract, contingency reserves, energy limits
      mdo.results.Pc(:,t)  = baseMVA * mdo.QP.x(vv.i1.Pc(t):vv.iN.Pc(t));
      mdo.results.Rpp(:,t) = baseMVA * mdo.QP.x(vv.i1.Rpp(t):vv.iN.Rpp(t));
      mdo.results.Rpm(:,t) = baseMVA * mdo.QP.x(vv.i1.Rpm(t):vv.iN.Rpm(t));
      if ns
        mdo.results.Sm(:,t)  = baseMVA * mdo.QP.x(vv.i1.Sm(t):vv.iN.Sm(t));
        mdo.results.Sp(:,t)  = baseMVA * mdo.QP.x(vv.i1.Sp(t):vv.iN.Sp(t));
      end
    end
    % Ramping reserves
    for t = 1:mdo.idx.ntramp
      mdo.results.Rrp(:,t) = baseMVA * mdo.QP.x(vv.i1.Rrp(t):vv.iN.Rrp(t));
      mdo.results.Rrm(:,t) = baseMVA * mdo.QP.x(vv.i1.Rrm(t):vv.iN.Rrm(t));
    end
    % Expected energy prices for generators, per generator and per period,
    % both absolute and conditional on making it to that period
    mdo.results.GenPrices = zeros(ng, nt);
    mdo.results.CondGenPrices = zeros(ng, nt);
    for t = 1:nt
      pp = zeros(ng,1);
      for j = 1:mdo.idx.nj(t)
        for k = 1:mdo.idx.nc(t,j)+1
          pp = pp + mdo.flow(t,j,k).mpc.bus(mdo.flow(t,j,k).mpc.gen(:,GEN_BUS), LAM_P);
        end
      end
      mdo.results.GenPrices(:,t) = pp;
      mdo.results.CondGenPrices(:, t) = pp / mdo.StepProb(t);
    end
    % Obtain contingency reserve prices, per generator and period
    mdo.results.RppPrices = zeros(ng, nt);
    mdo.results.RpmPrices = zeros(ng, nt);
    for t = 1:nt
      mdo.results.RppPrices(:, t) = mdo.QP.lambda.lower(vv.i1.Rpp(t):vv.iN.Rpp(t)) / baseMVA;
      mdo.results.RpmPrices(:, t) = mdo.QP.lambda.lower(vv.i1.Rpm(t):vv.iN.Rpm(t)) / baseMVA;
      for j = 1:mdi.idx.nj(t);
        for k = 1:mdi.idx.nc(t,j)+1
          ii = find(mdi.flow(t,j,k).mpc.gen(:, GEN_STATUS) > 0);
          mdo.results.RppPrices(ii, t) = mdo.results.RppPrices(ii, t) + mdo.QP.lambda.mu_l(ll.i1.dPpRp(t,j,k):ll.iN.dPpRp(t,j,k)) / baseMVA;
          mdo.results.RpmPrices(ii, t) = mdo.results.RpmPrices(ii, t) + mdo.QP.lambda.mu_l(ll.i1.dPmRm(t,j,k):ll.iN.dPmRm(t,j,k)) / baseMVA;
        end
      end
    end
    % Obtain ramping reserve prices, per generator and period
    mdo.results.RrpPrices = zeros(ng, mdo.idx.ntramp);
    mdo.results.RrmPrices = zeros(ng, mdo.idx.ntramp);
    % First, 1:nt-1
    for t = 1:nt-1
      for j1 = 1:mdo.idx.nj(t)
        for j2 = 1:mdo.idx.nj(t+1)
          if mdi.tstep(t+1).TransMask(j2,j1)
            mdo.results.RrpPrices(:, t) = mdo.results.RrpPrices(:, t) + mdo.QP.lambda.mu_l(ll.i1.Rrp(t,j1,j2):ll.iN.Rrp(t,j1,j2)) / baseMVA;
            mdo.results.RrmPrices(:, t) = mdo.results.RrmPrices(:, t) + mdo.QP.lambda.mu_l(ll.i1.Rrm(t,j1,j2):ll.iN.Rrm(t,j1,j2)) / baseMVA;
          end
        end
      end
    end
    % then last period only if specified for with terminal state
    if ~mdo.OpenEnded
      for j1 = 1:mdo.idx.nj(nt)
        mdo.results.RrpPrices(:, nt) = mdo.results.RrpPrices(:, nt) + mdo.QP.lambda.mu_l(ll.i1.Rrp(nt,j1,1):ll.iN.Rrp(nt,j1,1)) / baseMVA;
        mdo.results.RrmPrices(:, nt) = mdo.results.RrmPrices(:, nt) + mdo.QP.lambda.mu_l(ll.i1.Rrm(nt,j1,1):ll.iN.Rrm(nt,j1,1)) / baseMVA;
      end
    end
    % Expected wear and tear costs per gen and period
    mdo.results.ExpectedRampCost = zeros(ng, mdo.idx.ntramp+1);
    % First do first period wrt to InitialPg.
    for j = 1:mdi.idx.nj(1)
      w = mdo.tstep(1).TransMat(j,1); % the probability of going from initial state to jth
      mdo.results.ExpectedRampCost(:, 1) = mdo.results.ExpectedRampCost(:, 1) ...
          + 0.5 * w * mdo.RampWearCostCoeff(:,1) .* (mdo.flow(1,j,1).mpc.gen(:,PG) - mdo.InitialPg).^2;
    end
    % Then the remaining periods
    for t = 2:nt
      for j2 = 1:mdo.idx.nj(t)
        for j1 = 1:mdo.idx.nj(t-1)
          w = mdo.tstep(t).TransMat(j2,j1) * mdo.CostWeights(1, j1, t-1);
          mdo.results.ExpectedRampCost(:, t) = mdo.results.ExpectedRampCost(:, t) ...
              + 0.5 * w * mdo.RampWearCostCoeff(:,t) .* (mdo.flow(t,j2,1).mpc.gen(:,PG) - mdo.flow(t-1,j1,1).mpc.gen(:,PG)) .^2;
        end
      end
    end
    % Finally, if there is a terminal state problem, apply cost to
    if ~mdo.OpenEnded
      for j = 1:mdi.idx.nj(nt)
        w = mdi.tstep(t+1).TransMat(1, j) * mdi.CostWeights(1, j, nt);
        mdo.results.ExpectedRampCost(:, nt+1) = 0.5 * w * mdo.RampWearCostCoeff(:,nt+1) .* (mdo.TerminalPg - mdo.flow(nt,j,1).mpc.gen(:,PG)) .^2;
      end
    end
    % Compute expected dispatch, conditional on making it to the
    % corresponding period
    mdo.results.ExpectedDispatch = zeros(ng, nt);
    for t = 1:nt
      pp = sum(mdo.CostWeights(1,1:mdo.idx.nj(t),t)');    % gamma(t+1)
      for j = 1:mdo.idx.nj(t)
        mdo.results.ExpectedDispatch(:,t) = mdo.results.ExpectedDispatch(:,t) + ...
              mdo.CostWeights(1,j,t)/pp * mdo.flow(t,j,1).mpc.gen(:,PG);
      end
    end
    % If Cyclic storage, pull InitialStorage value out of x
    if ns && mdo.Storage.ForceCyclicStorage
      mdo.Storage.InitialStorage = baseMVA * mdo.QP.x(vv.i1.S0:vv.iN.S0);
    end
    % Compute expected storage state trajectory
    mdo.Storage.ExpectedStorageState = zeros(ns,nt);
    if ns
      for t = 1:nt
        pp = sum(mdo.CostWeights(1,1:mdo.idx.nj(t),t)');    %% gamma(t+1)
        for j = 1:mdo.idx.nj(t)
          Lfj = mdo.tstep(t).Lf( j:mdo.idx.nj(t):(ns-1)*mdo.idx.nj(t)+j, :);
          Ngj = mdo.tstep(t).Ng( j:mdo.idx.nj(t):(ns-1)*mdo.idx.nj(t)+j, :);
          Nhj = mdo.tstep(t).Nh( j:mdo.idx.nj(t):(ns-1)*mdo.idx.nj(t)+j, :);
          mdo.Storage.ExpectedStorageState(:,t) = ...
              mdo.Storage.ExpectedStorageState(:,t) + ...
                  baseMVA * mdo.CostWeights(1,j,t)/pp * ...
                  ( Lfj * mdo.Storage.InitialStorage/baseMVA + ...
                    (Ngj + Nhj) * mdo.QP.x );
        end
      end
      mdo.Storage.ExpectedStorageDispatch = ...
          mdo.results.ExpectedDispatch(mdo.Storage.UnitIdx, :);
    end
    % If there is a dynamical system, extract the state vectors and outputs
    % from the solution
    if ntds
      if nzds
        mdo.results.Z = zeros(nzds, ntds);
        for t = 1:ntds
          mdo.results.Z(:,t) = mdo.QP.x(vv.i1.Z(t):vv.iN.Z(t));
        end
      end
      mdo.results.Y = zeros(nyds, ntds);
      if nyds
        for t = 1:ntds
          mdo.results.Y(:, t) = ...
                  mdo.QP.A(ll.i1.DSy(t):ll.iN.DSy(t), :) * mdo.QP.x;
        end
      end
    end
    mdo.results.f = mdo.QP.f;
  end   % if success
end     % if mpopt.most.solve_model

mdo.results.success = success;
mdo.results.SetupTime = et_setup;
mdo.results.SolveTime = toc(t0);

if verbose
  fprintf('- MOST: Done.\n\n');
end

function Ostr = most(Istr, mpopt)
%MOST MATPOWER Optimal Scheduling Tool
%   MD_OUT = MOST(MD_IN)
%   MD_OUT = MOST(MD_IN, MPOPT)
%
%   Solves multi-period, stochastic, contingency constrained, optimal
%   power flow problem with linear constraints and unit commitment.
%   Depending on inputs it may include DC power flow constraints or
%   a simple total power balance condition.
%
%   Inputs:
%       MD_IN   MOST data structure, input
%           (see docs elsewhere for details)
%       MPOPT   MATPOWER options struct, relevant fields are (default
%               value in parens):
%           verbose - see 'help mpoption'
%           <solver specific options> - e.g. cplex, gurobi, etc,
%                     see 'help mpoption'
%           most.build_model (1) - build the MIQP, both constraints and
%                   standard costs (not coordination cost) and store in
%                   QP field of MD_OUT
%           most.solve_model (1) - solve the MIQP; if coordination
%                   cost exists, update it; requires either 'most.build_model'
%                   set to 1 or MD_IN.QP must contain previously built model
%           most.resolve_new_cost (0) - use when MIQP is already built and
%                   unchanged except for new coordination cost
%           most.dc_model (1) - use DC flow network model as opposed to simple
%                   generation = demand constraint
%           most.q_coordination (0) - create Qg variables for reactive power
%                   coordination
%           most.security_constraints (-1) - include contingency contstraints,
%                   if present, 1 = always include, 0 = never include
%           most.storage.terminal_target (-1) - constrain the expected terminal
%                   storage to target value, if present (1 = always, 0 = never)
%           most.storage.cyclic (0) - if 1, then initial storage is a variable
%                   constrained to = final expected storage; can't be
%                   simultaneously true with most.storage.terminal_target
%           most.uc.run (-1) - flag to indicate whether to perform unit
%                   commitment; 0 = do NOT perform UC, 1 = DO perform UC,
%                   -1 = perform UC if MD_IN.UC.CommitKey is present/non-empty
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

%   SuperOPF
%   $Id$
%   by Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Nacional de Colombia
%   and Ray Zimmerman, PSERC Cornell
%   Copyright (c) 2010-2015 by Power System Engineering Research Center (PSERC)

%   mpc is a MATPOWER case file or case struct with the fields baseMVA, bus,
%   gen, branch. It may also include a 'contingencies'
%   field (in place of contab argument) and the fields 'reserve' and
%   'energy_delta_cost' (in place of the offer argument). In this case, the
%   'reserve' and 'energy_delta_cost' fields look like the following, where
%   offerp refers to the first ng rows of offer and offerq to an optional
%   2nd set of ng rows:
%       .reserve
%           .cost
%               .Rp_pos     [ offerp(:, 1) ]
%               .Rp_neg     [ offerp(:, 3) ]
%               .Rq_pos     [ offerq(:, 1) ]    (optional)
%               .Rq_neg     [ offerq(:, 3) ]    (optional)
%           .cap
%               .Rp_pos     [ offerp(:, 2) ]
%               .Rp_neg     [ offerp(:, 4) ]
%               .Rq_pos     [ offerq(:, 2) ]    (optional)
%               .Rq_neg     [ offerq(:, 4) ]    (optional)
%       .energy_delta_cost
%           .dP_pos         [ offerp(:, 5) ]
%           .dP_neg         [ offerp(:, 6) ]
%           .dQ_pos         [ offerq(:, 5) ]    (optional)
%           .dQ_neg         [ offerq(:, 6) ]    (optional)
%
%   offer can also be a struct that contains the following fields:
%       offer.PositiveActiveReservePrice
%       offer.PositiveActiveReserveQuantity
%       offer.NegativeActiveReservePrice
%       offer.NegativeActiveReserveQuantity
%       offer.PositiveActiveDeltaPrice
%       offer.NegativeActiveDeltaPrice

% Schedule
%         .contingencies
%         .offer
%         .UC.CommitSched           ng x nt
%         .OpCondSched().tab
%         .Delta_T
%         .Initial.Pg0
%         .Initial.PgoPg1_Enforce
%         .Terminal_Pg
%         .Initial_Pg
%         .Storage.Type
%         .Storage.MinStorageLevel
%         .Storage.MaxStorageLevel
%         .Storage.InitialStorage
%         .Storage.OutEff
%         .Storage.InEff

% Istr
% To summarize, for (most.build_model, most.solve_model, most.resolve_new_cost) ...
%   (1,0,?) -> set up the problem, but don't solve it
%   (1,1,?) -> set up the problem, then solve it
%   (0,1,1) -> assume it's already set up, solve it with new costs
%   (0,1,0) -> replace with (1,1,0) internally
%   (0,0,?) -> error
% Istr.mpc                  initial system data baseMVA, bus, gen, branch
% Istr.UC.CommitSched       ng x nt
% Istr.Delta_T              length of time step in hours
% Istr.Augmented Lagrangian
% Istr.Storage.UnitIdx     Which gens are Storage Units
% Istr.Storage.Type        1) Pumped 2) Battery 3) Compressed air
% Istr.Storage.MinStorageLevel(t)   in MWH
% Istr.Storage.MaxStorageLevel(t)
% Istr.Storage.InitialStorage
% Istr.Storage.OutEff
% Istr.Storage.InEff
% Istr.idx.nt
% Istr.cont(t,j).contab
% Istr.tstep(t)            input data for each period
% Istr.tstep(t).offer
% Istr.tstep(t).OpCondSched(j).tab  Mods from Istr.mpc for each base scenario,
%                        except for generator commitment status which comes in Istr.UC.CommitSched
% Istr.tstep(t).TransMat(j2,j1)  Transition probabilities from j1(t-1) to j2(t)
% Istr.tstep(t).TransMask(j2,j1)  Ramp reserve mask for transition from j1(t-1) to j2(t)

% Istr.idx.nj(t)
% Istr.idx.nc(t,j)
% Istr.idx.nb_total
% Istr.idx.nb(t,j,k)
% Istr.idx.nf_total
% Istr.flow(t,j,k).mpc (with bus gen, branch, gencost)

tmptime(1,:) = clock;

%% default arguments
if nargin < 2
    mpopt = mpoption;       %% use default options
end

verbose = mpopt.verbose;

if verbose
    fprintf('\n=============================================================================\n');
    fprintf(  '          MATPOWER Optimal Scheduling Tool  --  MOST Version %s\n', mostver());
    fprintf(  '          A multi-period stochastic secure OPF with unit commitment\n');
    fprintf(  '                       -----  Built on MATPOWER  -----\n');
    fprintf(  '  by Carlos E. Murillo-Sanchez, Universidad Nacional de Colombia--Manizales\n');
    fprintf(  '                  and Ray D. Zimmerman, Cornell University\n');
    fprintf(  '                 (c) 2012-2016 Cornell University and PSERC\n');
    fprintf(  '=============================================================================\n');
end

%% if you want to do a normal solve, you have to create the QP
if mpopt.most.solve_model && ~mpopt.most.resolve_new_cost
  mpopt = mpoption(mpopt, 'most.build_model', 1);
end
if ~mpopt.most.build_model && ~mpopt.most.solve_model
  error('most: Ah ... are you sure you want to do nothing? (either ''most.build_model'' or ''most.solve_model'' must be true)');
end
if ~isfield(Istr, 'IncludeFixedReserves')
  Istr.IncludeFixedReserves = false;
end

%% set up some variables we use throughout
ng = size(Istr.mpc.gen, 1);
nt = Istr.idx.nt;
ns = length(Istr.Storage.UnitIdx);
Istr.idx.ng = ng;
Istr.idx.ns = ns;
baseMVA = Istr.mpc.baseMVA;
for t = 1:nt
  Istr.idx.nj(t) = length(Istr.tstep(t).OpCondSched);
end
if ~isfield(Istr, 'UC') || ~isfield(Istr.UC, 'CommitSched') || ...
        isempty(Istr.UC.CommitSched)
  if isfield(Istr, 'CommitSched') && ~isempty(Istr.CommitSched)
    warning('-----  most: MD_IN.CommitSched has moved to MD_IN.UC.CommitSched, please update your code.  -----');
    Istr.UC.CommitSched = Istr.CommitSched;
  else
    error('most: commitment schedule must be provided in MSPD_IN.UC.CommitSched');
  end
end

% set up model options
UC = mpopt.most.uc.run;
if UC == -1
  if isempty(Istr.UC.CommitKey)
    UC = 0;
  else
    UC = 1;
  end
end
mo = struct(...
    'DCMODEL',                      mpopt.most.dc_model, ...
    'SecurityConstrained',          mpopt.most.security_constraints, ...
    'QCoordination',                mpopt.most.q_coordination, ...
    'ForceCyclicStorage',           mpopt.most.storage.cyclic, ...
    'CyclicCommitment',             mpopt.most.uc.cyclic, ...
    'ForceExpectedTerminalStorage', mpopt.most.storage.terminal_target, ...
    'alpha',                        mpopt.most.alpha ...
);
if mo.SecurityConstrained == -1
  if isfield(Istr, 'cont') && isfield(Istr.cont(1,1), 'contab') && ...
          ~isempty(Istr.cont(1,1).contab)
    mo.SecurityConstrained = 1;
  else
    mo.SecurityConstrained = 0;
  end
end
if mo.ForceExpectedTerminalStorage == -1
  if isempty(Istr.Storage.ExpectedTerminalStorageAim) && ...
     isempty(Istr.Storage.ExpectedTerminalStorageMin) && ...
     isempty(Istr.Storage.ExpectedTerminalStorageMax)
    mo.ForceExpectedTerminalStorage = 0;
  else
    mo.ForceExpectedTerminalStorage = 1;
  end
end
if mo.SecurityConstrained && ~(isfield(Istr, 'cont') && ...
      isfield(Istr.cont(1,1), 'contab') && ~isempty(Istr.cont(1,1).contab))
  error('most: MD_IN.cont(t,j).contab cannot be empty when MPOPT.most.security_constraints = 1');
end
if mo.ForceExpectedTerminalStorage == 1;
  if mo.ForceCyclicStorage
    error('most: storage model cannot be both cyclic and include a terminal target value; must change MPOPT.most.storage.cyclic or MPOPT.most.storage.terminal_target');
  end
  if ns && isempty(Istr.Storage.ExpectedTerminalStorageAim) && ...
           isempty(Istr.Storage.ExpectedTerminalStorageMin) && ...
           isempty(Istr.Storage.ExpectedTerminalStorageMax)
    error('most: MD_IN.Storage.ExpectedTerminalStorageAim|Min|Max cannot all be empty when MPOPT.most.storage.terminal_target = 1');
  end
  if ~isempty(Istr.Storage.ExpectedTerminalStorageAim)
    Istr.Storage.ExpectedTerminalStorageMin = Istr.Storage.ExpectedTerminalStorageAim;
    Istr.Storage.ExpectedTerminalStorageMax = Istr.Storage.ExpectedTerminalStorageAim;
  end
end
if UC && (~isfield(Istr.UC, 'CommitKey') || isempty(Istr.UC.CommitKey))
  error('most: cannot run unit commitment without specifying MD.UC.CommitKey');
end

if ns
  if isempty(Istr.Storage.InitialStorage)
    error('most: Storage.InitialStorage must be specified');
  end
  if isempty(Istr.Storage.TerminalChargingPrice0)
    Istr.Storage.TerminalChargingPrice0 = Istr.Storage.TerminalStoragePrice;
  end
  if isempty(Istr.Storage.TerminalDischargingPrice0)
    Istr.Storage.TerminalDischargingPrice0 = Istr.Storage.TerminalStoragePrice;
  end
  if isempty(Istr.Storage.TerminalChargingPriceK)
    Istr.Storage.TerminalChargingPriceK = Istr.Storage.TerminalStoragePrice;
  end
  if isempty(Istr.Storage.TerminalDischargingPriceK)
    Istr.Storage.TerminalDischargingPriceK = Istr.Storage.TerminalStoragePrice;
  end
  if isempty(Istr.Storage.MinStorageLevel)
    error('most: Storage.MinStorageLevel must be specified');
  else
    MinStorageLevel = Istr.Storage.MinStorageLevel;
  end
  if size(MinStorageLevel, 1) == 1 && ns > 1    %% expand rows
    MinStorageLevel = ones(ns, 1) * MinStorageLevel;
  end
  if size(MinStorageLevel, 2) == 1 && nt > 1    %% expand cols
    MinStorageLevel = MinStorageLevel * ones(1, nt);
  end
  if isempty(Istr.Storage.MaxStorageLevel)
    error('most: Storage.MaxStorageLevel must be specified');
  else
    MaxStorageLevel = Istr.Storage.MaxStorageLevel;
  end
  if size(MaxStorageLevel, 1) == 1 && ns > 1    %% expand rows
    MaxStorageLevel = ones(ns, 1) * MaxStorageLevel;
  end
  if size(MaxStorageLevel, 2) == 1 && nt > 1    %% expand cols
    MaxStorageLevel = MaxStorageLevel * ones(1, nt);
  end
  if isempty(Istr.Storage.InEff)
    InEff = 1;                      %% no efficiency loss by default
  else
    InEff = Istr.Storage.InEff;
  end
  if size(InEff, 1) == 1 && ns > 1  %% expand rows
    InEff = ones(ns, 1) * InEff;
  end
  if size(InEff, 2) == 1 && nt > 1  %% expand cols
    InEff = InEff * ones(1, nt);
  end
  if isempty(Istr.Storage.OutEff)
    OutEff = 1;                     %% no efficiency loss by default
  else
    OutEff = Istr.Storage.OutEff;
  end
  if size(OutEff, 1) == 1 && ns > 1 %% expand rows
    OutEff = ones(ns, 1) * OutEff;
  end
  if size(OutEff, 2) == 1 && nt > 1 %% expand cols
    OutEff = OutEff * ones(1, nt);
  end
  if isempty(Istr.Storage.LossFactor)
    LossFactor = 0;                     %% no losses by default
  else
    LossFactor = Istr.Storage.LossFactor;
  end
  if size(LossFactor, 1) == 1 && ns > 1 %% expand rows
    LossFactor = ones(ns, 1) * LossFactor;
  end
  if size(LossFactor, 2) == 1 && nt > 1 %% expand cols
    LossFactor = LossFactor * ones(1, nt);
  end
  if isempty(Istr.Storage.rho)
    rho = 1;                        %% use worst case by default (for backward compatibility)
  else
    rho = Istr.Storage.rho;
  end
  if size(rho, 1) == 1 && ns > 1    %% expand rows
    rho = ones(ns, 1) * rho;
  end
  if size(rho, 2) == 1 && nt > 1    %% expand cols
    rho = rho * ones(1, nt);
  end
  if isempty(Istr.Storage.InitialStorageLowerBound)
    if mo.ForceCyclicStorage        %% lower bound for var s0, take from t=1
      Istr.Storage.InitialStorageLowerBound = MinStorageLevel(:, 1);
    else                            %% Sm(0), default = fixed param s0
      Istr.Storage.InitialStorageLowerBound = Istr.Storage.InitialStorage;
    end
  end
  if isempty(Istr.Storage.InitialStorageUpperBound)
    if mo.ForceCyclicStorage        %% upper bound for var s0, take from t=1
      Istr.Storage.InitialStorageUpperBound = MaxStorageLevel(:, 1);
    else                            %% Sp(0), default = fixed param s0
      Istr.Storage.InitialStorageUpperBound = Istr.Storage.InitialStorage;
    end
  end
  
  LossCoeff = Istr.Delta_T * LossFactor/2;
  beta1 = (1-LossCoeff) ./ (1+LossCoeff);
  beta2 = 1 ./ (1+LossCoeff);
  beta3 = 1 ./ (1/(1-mo.alpha) + LossCoeff);
  beta4 = mo.alpha/(1-mo.alpha) * beta2 .* beta3;
  beta5 = beta1 ./ beta2 .* (beta3 + beta4);
  beta2EtaIn      = Istr.Delta_T * beta2 .* InEff;
  beta2overEtaOut = Istr.Delta_T * beta2 ./ OutEff;
  beta3EtaIn      = Istr.Delta_T * beta3 .* InEff;
  beta3overEtaOut = Istr.Delta_T * beta3 ./ OutEff;
  beta4EtaIn      = Istr.Delta_T * beta4 .* InEff;
  beta4overEtaOut = Istr.Delta_T * beta4 ./ OutEff;
  diagBeta2EtaIn1      = spdiags(beta2EtaIn(:,1),      0, ns, ns);
  diagBeta2overEtaOut1 = spdiags(beta2overEtaOut(:,1), 0, ns, ns);
  diagBeta3EtaIn1      = spdiags(beta3EtaIn(:,1),      0, ns, ns);
  diagBeta3overEtaOut1 = spdiags(beta3overEtaOut(:,1), 0, ns, ns);
  diagBeta4EtaIn1      = spdiags(beta4EtaIn(:,1),      0, ns, ns);
  diagBeta4overEtaOut1 = spdiags(beta4overEtaOut(:,1), 0, ns, ns);
end
if ~isfield(Istr.idx, 'nyt') || isempty(Istr.idx.nyt) || ~Istr.idx.nyt
  nyt = 0;
  nzd = 0;
  nyo = 0;
else
  nyt = Istr.idx.nyt;
  nzd = size(Istr.dstep(1).A, 1);
  nyo = size(Istr.dstep(1).D, 1);   % # of outputs of dynamical system
end
Istr.idx.nyt = nyt;
Istr.idx.nzd = nzd;
Istr.idx.nyo = nyo;

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

% Make data tables with full # of cols and add also pseudo OPF results to
% be able to run printpf on them
Istr.mpc.bus(:, MU_VMIN) = 0;
Istr.mpc.gen(:, MU_QMIN) = 0;
Istr.mpc.branch(:,MU_ANGMAX) = 0;
Istr.mpc.f = 0;
Istr.mpc.et = 0;
Istr.mpc.success = 1;

if mpopt.most.build_model
  if verbose
    fprintf('- Building indexing structures.\n');
  end

  %% save model options in data structure
  Istr.DCMODEL                              = mo.DCMODEL;
  Istr.SecurityConstrained                  = mo.SecurityConstrained;
  Istr.QCoordination                        = mo.QCoordination;
  Istr.Storage.ForceCyclicStorage           = mo.ForceCyclicStorage;
  Istr.Storage.ForceExpectedTerminalStorage = mo.ForceExpectedTerminalStorage;
  Istr.UC.run                               = UC;
  Istr.UC.CyclicCommitment                  = mo.CyclicCommitment;
  Istr.alpha                                = mo.alpha;
  if ~isfield(Istr, 'OpenEnded'), Istr.OpenEnded = 1; end

  if UC
    % Make sure MinUp and MinDown are all >= 1
    if any(Istr.UC.MinUp < 1) && any(Istr.UC.MinUp < 1)
        error('most: UC.MinUp and UC.MinDown must all be >= 1');
    end
    % Unless something is forced off in Istr.CommitKey, or as a result of
    % not fulfilling its Istr.UC.MinDown in early periods, it should be available
    % for commitment and thus a contingency including its outage should not
    % be deleted.
    Istr.UC.CommitSched = (Istr.UC.CommitKey >= 0);  % Treat anything but -1 as on.
    if ~Istr.UC.CyclicCommitment
      for i = 1:ng
        if Istr.UC.InitialState(i) < 0
          nn = Istr.UC.MinDown(i) + Istr.UC.InitialState(i); % time to go before startup
          if nn > 0
            Istr.UC.CommitSched(i, 1:nn) = 0;
          end
        elseif Istr.UC.InitialState(i) > 0
          nn = Istr.UC.MinUp(i) - Istr.UC.InitialState(i); % time to go before shutdown
          if nn > 0
            Istr.UC.CommitSched(i, 1:nn) = 1;
          end
        end
      end
    end
  end
  % From now on, Istr.UC.CommitSched has zeros for stuff that is definitely
  % off, and ones for stuff that might be on, so those zeros can be used to
  % trim off contingencies that won't happen.
  % Start by creating the base flow data for all scenarios
  for t = 1:nt
    for j = 1:Istr.idx.nj(t)
      mpc = Istr.mpc;
      mpc.gen(:, GEN_STATUS) = Istr.UC.CommitSched(:, t);

      %% for backward compatibility, putting time dependent energy offer
      %% data in offer(t).gencost is deprecated, please use profiles
      if isfield(Istr.offer(t), 'gencost') && ~isempty(Istr.offer(t).gencost)
        mpc.gencost = Istr.offer(t).gencost;
      end

      if ~isempty(Istr.tstep(t).OpCondSched(j).tab)
        changelist = unique(Istr.tstep(t).OpCondSched(j).tab(:, CT_LABEL));
        for label = changelist'
          mpc = apply_changes(label, mpc, Istr.tstep(t).OpCondSched(j).tab);
        end
      end
      Istr.flow(t,j,1).mpc = mpc;
      Istr.idx.nb(t,j,1) = size(Istr.flow(t,j,1).mpc.bus, 1);
      Istr.idx.ny(t,j,1) = length(find(Istr.flow(t,j,1).mpc.gencost(:, MODEL) == PW_LINEAR));
    end
  end
  % Then continue to create contingent flow scenarios, deleting any
  % redundant contingencies (i.e., decommitting a gen or branch when its
  % status is guaranteed to be off). No rows are deleted from gen or branch,
  % but the number of contingencies can indeed change.
  for t = 1:nt
    % Set default ramp reserve mask, if not provided
    if ~isfield(Istr.tstep(t), 'TransMask') || isempty(Istr.tstep(t).TransMask)
      Istr.tstep(t).TransMask = ones(size(Istr.tstep(t).TransMat));
    end
    % First get current step's scenario probabilities
    if t == 1
      scenario_probs = Istr.tstep(1).TransMat; % the probability of the initial state is 1
    else
      scenario_probs = Istr.tstep(t).TransMat * Istr.CostWeights(1, 1:Istr.idx.nj(t-1), t-1)'; % otherwise compute from previous step base cases
    end
    Istr.StepProb(t) = sum(scenario_probs); % probability of making it to the t-th step
    if Istr.SecurityConstrained
      for j = 1:Istr.idx.nj(t)
        [tmp, ii] = sort(Istr.cont(t,j).contab(:, CT_LABEL)); %sort in ascending contingency label
        contab = Istr.cont(t,j).contab(ii, :);
        rowdecomlist = ones(size(contab,1), 1);
        for l = 1:size(contab, 1)
          if contab(l, CT_TABLE) == CT_TGEN  && contab(l, CT_COL) == GEN_STATUS ...
              && contab(l, CT_CHGTYPE) == CT_REP && contab(l, CT_NEWVAL) == 0 ... % gen turned off
              && Istr.flow(t,j,1).mpc.gen(contab(l, CT_ROW), GEN_STATUS) <= 0    % but it was off on input
           rowdecomlist(l) = 0;
          elseif contab(l, CT_TABLE) == CT_TBRCH && contab(l, CT_COL) == BR_STATUS ...
              && contab(l, CT_CHGTYPE) == CT_REP && contab(l, CT_NEWVAL) == 0 ... % branch taken out
              && Istr.flow(t,j,1).mpc.branch(contab(l, CT_ROW), BR_STATUS) <= 0  % but it was off on input
            rowdecomlist(l) = 0;
          end
        end
        contab = contab(rowdecomlist ~= 0, :);
        Istr.cont(t, j).contab = contab;
        clist = unique(contab(:, CT_LABEL));
        Istr.idx.nc(t, j) = length(clist);
        k = 2;
        for label = clist'
          Istr.flow(t, j, k).mpc = apply_changes(label, Istr.flow(t, j, 1).mpc, contab);
          ii = find( label == contab(:, CT_LABEL) );
          Istr.CostWeights(k, j, t) = contab(ii(1), CT_PROB);
          Istr.idx.nb(t, j, k) = size(Istr.flow(t, j, k).mpc.bus, 1);
          Istr.idx.ny(t, j, k) = length(find(Istr.flow(t, j, 1).mpc.gencost(:, MODEL) == PW_LINEAR));
          k = k + 1;
        end
        Istr.CostWeights(1, j, t) = 1 - sum(Istr.CostWeights(2:Istr.idx.nc(t,j)+1, j, t));
        Istr.CostWeights(1:Istr.idx.nc(t,j)+1, j, t) = scenario_probs(j) * Istr.CostWeights(1:Istr.idx.nc(t,j)+1, j, t);
      end
    else
      for j = 1:Istr.idx.nj(t)
        Istr.idx.nc(t, j) = 0;
        Istr.CostWeights(1, j, t) = scenario_probs(j);
      end
    end
  end

  % Compute adjusted (for alpha) cost weights for objective function
  if Istr.SecurityConstrained && Istr.alpha ~= 0
    for t = 1:nt
      for j = 1:Istr.idx.nj(t)
        Istr.CostWeightsAdj(1, j, t) = Istr.CostWeights(1, j, t);
        for k = 2:Istr.idx.nc(t,j)+1
          Istr.CostWeightsAdj(k, j, t) = (1-Istr.alpha) * Istr.CostWeights(k, j, t);
          Istr.CostWeightsAdj(1, j, t) = Istr.CostWeightsAdj(1, j, t) + Istr.alpha * Istr.CostWeights(k, j, t);
        end
      end
    end
  else
    Istr.CostWeightsAdj = Istr.CostWeights;
  end

  % If UC, also need to (possibly) modify gencosts so that each fm(p) at
  % p=0 is zero, so that fm(p) + u*c00 is equal to the original f(p). This
  % portion of the cost is valid at t if the unit is commited there, but
  % only for base scenarios and contingencies in which this particular unit
  % is not ousted, so must be careful later when probability-weighting the
  % corresponding u(i,t)!
  if UC
    if ~isfield(Istr.UC, 'c00') || isempty(Istr.UC.c00) % if not empty assume
      Istr.UC.c00 = zeros(ng, nt);                      % contains correct info!
      for t = 1:nt
        for j = 1:Istr.idx.nj(t)
          for k = 1:Istr.idx.nc(t,j)+1
            mpc = Istr.flow(t,j,k).mpc;
            c00tjk = totcost(mpc.gencost, zeros(ng,1));
            Istr.UC.c00(:, t) = Istr.UC.c00(:, t) + Istr.CostWeightsAdj(k, j, t) * c00tjk;
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
            Istr.flow(t,j,k).mpc = mpc;
          end
        end
      end
    end
  end
  
  % Build variable indexing mechanism
  % Find total number of flows, buses and ny variables; (including offline gens)
  Istr.idx.nf_total = 0;
  Istr.idx.nb_total = 0;
  for t = 1:nt
    for j = 1:Istr.idx.nj(t)
      Istr.idx.nf_total = Istr.idx.nf_total + (1 + Istr.idx.nc(t,j));
      for k = 1:Istr.idx.nc(t,j)+1
        Istr.idx.nb_total = Istr.idx.nb_total + size(Istr.flow(t, j, k).mpc.bus, 1);
        ii = find(Istr.flow(t,j,k).mpc.gencost(:, MODEL) == PW_LINEAR);
        Istr.idx.ny(t,j,k) = length(ii);
      end
    end
  end
  Istr.idx.ns_total = ns * Istr.idx.nf_total;
  % Variable order resembles that of several C3SOPFs stacked together,
  % including the internally generated y variables, and then all of the
  % other new variables that are specific to HP, but excluding qg on one hand, and pc,
  % rp and rm since these now are common across several scenarios. So create first a matrix
  % of indices to the beginning of each c3sopf cell's vars.  Include the
  % mechanism for adding theta variables if we want to create DC flow restrictions.
  % Then start assigning the start and end indices for variables in each
  % c3sopf cell
  om = opt_model;
  nj_max = max(Istr.idx.nj);
  nc_max = max(max(Istr.idx.nc));
  Ing = speye(ng);
  Ins = speye(ns);
  if Istr.DCMODEL
    om = add_vars(om, 'Va', {nt, nj_max, nc_max+1});
  end
  om = add_vars(om, 'Pg', {nt, nj_max, nc_max+1});
  om = add_vars(om, 'dPp', {nt, nj_max, nc_max+1});
  om = add_vars(om, 'dPm', {nt, nj_max, nc_max+1});
  om = add_vars(om, 'y', {nt, nj_max, nc_max+1});
  if Istr.IncludeFixedReserves
    om = add_vars(om, 'R', {nt, nj_max, nc_max+1});
  end
  for t = 1:nt
    for j = 1:Istr.idx.nj(t)
      % first all angles if using DCMODEL
      for k = 1:Istr.idx.nc(t,j)+1
        if Istr.DCMODEL
          iref = find(Istr.flow(t,j,k).mpc.bus(:,BUS_TYPE) == REF);
          if verbose && length(iref) > 1
            errstr = ['\nmpsopfl: Warning: Multiple reference buses.\n', ...
              '           For a system with islands, a reference bus in each island\n', ...
              '           may help convergence, but in a fully connected system such\n', ...
              '           a situation is probably not reasonable.\n\n' ];
            fprintf(errstr);
          end
          Va0 = Istr.flow(t,j,k).mpc.bus(:,VA)*pi/180;
          Va_max = Inf(Istr.idx.nb(t,j,k), 1);
          Va_min = -Va_max;
          Va_min(iref) = Istr.flow(t,j,k).mpc.bus(iref,VA)*pi/180;
          Va_max(iref) = Va_min(iref);
          
          om = add_vars(om, 'Va', {t,j,k}, Istr.idx.nb(t,j,k), Va0, Va_min, Va_max);
        end
      end
      % All active injections in c3sopf cell
      for k = 1:Istr.idx.nc(t,j)+1
        mpc = Istr.flow(t,j,k).mpc;
        genmask = mpc.gen(:,GEN_STATUS) > 0;
        p0 = genmask .* mpc.gen(:,PG) / baseMVA;
        if UC       % relax bounds here, enforced by uPmax, uPmin constraints
          pmin = genmask .* (min(mpc.gen(:, PMIN) / baseMVA, 0) - 1);
          pmax = genmask .* (max(mpc.gen(:, PMAX) / baseMVA, 0) + 1);
        else        % enforce bounds here, subject to flow's GEN_STATUS
          pmin = genmask .* mpc.gen(:, PMIN) / baseMVA;
          pmax = genmask .* mpc.gen(:, PMAX) / baseMVA;
        end
        om = add_vars(om, 'Pg', {t,j,k}, ng, p0, pmin, pmax);
      end
      if Istr.IncludeFixedReserves
        for k = 1:Istr.idx.nc(t,j)+1
          mpc = Istr.flow(t,j,k).mpc;
          r = Istr.FixedReserves(t,j,k);
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
          Istr.FixedReserves(t,j,k).rgens = r.rgens;
          Istr.FixedReserves(t,j,k).igr   = r.igr;
          if isfield(r, 'original')
              Istr.FixedReserves(t,j,k).original = r.original;
          end
          Istr.FixedReserves(t,j,k)       = r;  %% for cost & qty (now that fields match)
          Rmax = Inf(ngr, 1);               %% bound above by ...
          kk = find(mpc.gen(r.igr, RAMP_10));
          Rmax(kk) = mpc.gen(r.igr(kk), RAMP_10);   %% ... ramp rate and ...
          kk = find(r.qty(r.igr) < Rmax);
          Rmax(kk) = r.qty(r.igr(kk));      %% ... stated max reserve qty
          Rmax = Rmax / baseMVA;
          om = add_vars(om, 'R', {t,j,k}, ngr, [], zeros(ngr, 1), Rmax);
        end
      end
      % All deltaP plus in c3sopf cell
      for k = 1:Istr.idx.nc(t,j)+1
        om = add_vars(om, 'dPp', {t,j,k}, ng, [], zeros(ng,1), []);
      end
      % All deltaP minus in c3sopf cell
      for k = 1:Istr.idx.nc(t,j)+1
        om = add_vars(om, 'dPm', {t,j,k}, ng, [], zeros(ng,1), []);
      end
      % All y variables in c3sopf cell - even if not committed.  There must
      % be a fixed cost associated with u(t,i,j) such that if u(t,i,j) = 0,
      % then the cost interpolated from the (x,y) pairs is zero, and if
      % u(t,i,j) = 1, then the fixed cost plus that interpolated from the
      % (x,y) pairs is as desired.
      for k = 1:Istr.idx.nc(t,j)+1
        om = add_vars(om, 'y', {t,j,k}, Istr.idx.ny(t,j,k), [], [], []);
      end %
    end % for j
  end % for t
  % Continue with pc, rpp, rpm, one set for each time period
  om = add_vars(om, 'Pc', {nt});
  om = add_vars(om, 'Rpp', {nt});
  om = add_vars(om, 'Rpm', {nt});
  for t = 1:nt
    om = add_vars(om, 'Pc', {t}, ng);
    %% non-negativity on Rpp and Rpm is redundant, leave unbounded below
    %% (except where gen is off-line for all j and k)
    Rpmin = -Inf(ng,1);
    off = ones(ng,1);
    for j = 1:Istr.idx.nj(t);
      for k = 1:Istr.idx.nc(t,j)+1
        off = off & Istr.flow(t,j,k).mpc.gen(:, GEN_STATUS) <= 0;
      end
    end
    Rpmin(off == 1) = 0;
    om = add_vars(om, 'Rpp', {t}, ng, [], Rpmin, Istr.offer(t).PositiveActiveReserveQuantity/baseMVA);
    om = add_vars(om, 'Rpm', {t}, ng, [], Rpmin, Istr.offer(t).NegativeActiveReserveQuantity/baseMVA);
  end
  % Now load following ramping reserves.  In open ended problem, we need to
  % specify nt-1 ramping reserves, those needed to transition 1-2, 2-3, ..
  % (nt-1)-nt . The initial ramp constraint (from t=0 to t=1) is data, not a
  % variable.  But in terminal state (at t=nt+1) or cyclical problems (when
  % the t=nt to t=1 transition is also considered) we need nt ramping
  % reserves.
  if ~Istr.OpenEnded
    Istr.idx.ntramp = nt;
  else
    Istr.idx.ntramp = nt - 1;
  end
  om = add_vars(om, 'Rrp', {Istr.idx.ntramp});
  om = add_vars(om, 'Rrm', {Istr.idx.ntramp});
  for t = 1:Istr.idx.ntramp
    ramp30 = Istr.flow(t,1,1).mpc.gen(:,RAMP_30)*2*Istr.Delta_T;
    om = add_vars(om, 'Rrp', {t}, ng, [], zeros(ng,1), ...
        min(Istr.offer(t).PositiveLoadFollowReserveQuantity, ramp30)/baseMVA);
%     om = add_vars(om, 'Rrm', {t}, ng, [], zeros(ng,1), ...
%         min(Istr.offer(t).NegativeLoadFollowReserveQuantity, ramp30)/baseMVA);
  end
  for t = 1:Istr.idx.ntramp
    ramp30 = Istr.flow(t,1,1).mpc.gen(:,RAMP_30)*2*Istr.Delta_T;
    om = add_vars(om, 'Rrm', {t}, ng, [], zeros(ng,1), ...
        min(Istr.offer(t).NegativeLoadFollowReserveQuantity, ramp30)/baseMVA);
  end
  % Continue with storage charge/discharge injections, one of each
  % for each flow; first all charge injections, then all discharge
  % injections
  om = add_vars(om, 'Psc', {nt, nj_max, nc_max+1});
  om = add_vars(om, 'Psd', {nt, nj_max, nc_max+1});
  for t = 1:nt
    for j = 1:Istr.idx.nj(t)
      for k = 1:Istr.idx.nc(t,j)+1
        if ns
          om = add_vars(om, 'Psc', {t,j,k}, ns, [], [], zeros(ns,1));
%           om = add_vars(om, 'Psd', {t,j,k}, ns, [], zeros(ns,1), []);
        end
      end
    end
  end
  if ns
    for t = 1:nt
      for j = 1:Istr.idx.nj(t)
        for k = 1:Istr.idx.nc(t,j)+1
          om = add_vars(om, 'Psd', {t,j,k}, ns, [], zeros(ns,1), []);
        end
      end
    end
  end
  % Continue with storage upper and lower bounds, one for each time period
  % and unit
  om = add_vars(om, 'Sp', {nt});
  om = add_vars(om, 'Sm', {nt});
  if ns
    for t = 1:nt
      om = add_vars(om, 'Sp', {t}, ns, [], [], MaxStorageLevel(:,t)/baseMVA);
%       om = add_vars(om, 'Sm', {t}, ns, [], MinStorageLevel(:,t)/baseMVA, []);
    end
    for t = 1:nt
      om = add_vars(om, 'Sm', {t}, ns, [], MinStorageLevel(:,t)/baseMVA, []);
    end
  end
  % Possible initial storage quantities when using cyclic storage dispatch
  % so that initial storage = expected terminal storage is a constraint
  if ns && Istr.Storage.ForceCyclicStorage
    om = add_vars(om, 'S0', ns, [], ...
        Istr.Storage.InitialStorageLowerBound / baseMVA, ...
        Istr.Storage.InitialStorageUpperBound / baseMVA);
  end
  % If there is a dynamical system with non-null state vector,
  % add those states here
  if nzd
    om = add_vars(om, 'Z', {nyt});
    for t = 1:nyt
      if t == 1
        zmin = Istr.z1;
        zmax = Istr.z1;
      else
        zmin = Istr.dstep(t).zmin;
        zmax = Istr.dstep(t).zmax;
      end
      z0 = (zmax - zmin) / 2;
      om = add_vars(om, 'Z', {t}, nzd, z0, zmin, zmax);
    end
  end
  % Now the integer variables; u variables mean on/off status
  if UC
    om = add_vars(om, 'u', {nt});
    om = add_vars(om, 'v', {nt});
    om = add_vars(om, 'w', {nt});
    vt0 = char('B' * ones(1, ng));  % default variable type for u is binary
    for t = 1:nt
      umin = zeros(ng, 1);
      umax = ones(ng, 1);
      % min up/down restrictions on u
      if ~Istr.UC.CyclicCommitment
        % min up time has not passed yet since startup occured, force ON
        umin( (Istr.UC.InitialState > 0) & ...
              (t+Istr.UC.InitialState-Istr.UC.MinUp <= 0) ) = 1;
        % min down time has not passed yet since shutdown occured, force OFF
        umax( (Istr.UC.InitialState < 0) & ...
              (t-Istr.UC.InitialState-Istr.UC.MinDown <= 0) ) = 0;
      end
      % set limits for units forced ON or forced OFF
      iON = find(Istr.UC.CommitKey(:,t) == 2);
      iOFF = find(Istr.UC.CommitKey(:,t) == -1);
      umin(iON)  = 1;
      umax(iOFF) = 0;

      % set variable types
      vt = vt0;                 % initialize all variable types to binary
      vt(umin == umax) = 'C';   % make continuous for those that are fixed
      
      om = add_vars(om, 'u', {t}, ng, zeros(ng, 1), umin, umax, vt);
    end
    % v variables mean startup events
    for t = 1:nt
      om = add_vars(om, 'v', {t}, ng, zeros(ng, 1), zeros(ng, 1), ones(ng, 1));
    end
    % w variables mean shutdown events
    for t = 1:nt
      om = add_vars(om, 'w', {t}, ng, zeros(ng, 1), zeros(ng, 1), ones(ng, 1));
    end
  end
  % An external program may be using coordination with AC flows, and in
  % that case we need corresponding Qg variables, whose only function is to
  % be constrained to zero if the commitment decision asks for that.
  if Istr.QCoordination
    om = add_vars(om, 'Qg', {nt, nj_max, nc_max+1});
    for t = 1:nt
      for j = 1:Istr.idx.nj(t)
        for k = 1:Istr.idx.nc(t,j)+1
          mpc = Istr.flow(t,j,k).mpc;
          genmask = mpc.gen(:,GEN_STATUS) > 0;
          q0 = genmask .* mpc.gen(:, QG) / baseMVA;
          if UC     % relax bounds here, enforced by uQmax, uQmin constraints
            qmin = genmask .* (min(mpc.gen(:, QMIN) / baseMVA, 0) - 1);
            qmax = genmask .* (max(mpc.gen(:, QMAX) / baseMVA, 0) + 1);
          else      % enforce bounds here, subject to flow's GEN_STATUS
            qmin = genmask .* mpc.gen(:, QMIN) / baseMVA;
            qmax = genmask .* mpc.gen(:, QMAX) / baseMVA;
          end
          om = add_vars(om, 'Qg', {t,j,k}, ng, q0, qmin, qmax);
        end
      end
    end
  end
  nvars = getN(om, 'var');
  Istr.idx.nvars = nvars;

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
    vv = get_idx(om);
    for t = 1:nt
      nsxnjt = ns*Istr.idx.nj(t);
      % Form G(t), H(t), B1(t), B2(t)
      G = sparse(nsxnjt, nvars);
      H = sparse(nsxnjt, nvars);
      B1 = sparse(nsxnjt, nsxnjt);
      B2 = sparse(nsxnjt, nsxnjt);
      for j = 1:Istr.idx.nj(t)
        ii  = ((1:ns)'-1)*Istr.idx.nj(t)+j;
        jj1 = (vv.i1.Psc(t,j,1):vv.iN.Psc(t,j,1))';
        jj2 = (vv.i1.Psd(t,j,1):vv.iN.Psd(t,j,1))';
        G = G + sparse(ii, jj1, -Istr.Delta_T  *  InEff(:,t), nsxnjt, nvars);
        H = H + sparse(ii, jj2, -Istr.Delta_T ./ OutEff(:,t), nsxnjt, nvars);
        B1 = B1 + sparse(ii, ii, beta1(:,t), nsxnjt, nsxnjt);
        B2 = B2 + sparse(ii, ii, beta2(:,t), nsxnjt, nsxnjt);
      end
      if t == 1
        % form Li, Lf, Mg, Mh, Ng, Nh, B1, B2 for t == 1
        jlist = [];
        for i=1:ns
          jlist = [ jlist; i*ones(Istr.idx.nj(t),1) ];
        end
        Istr.tstep(t).Li  = sparse((1:nsxnjt)', jlist, 1, nsxnjt, ns);
        Istr.tstep(t).Lf  = B1 * Istr.tstep(t).Li;
        Istr.tstep(t).Mg  = sparse(nsxnjt, nvars);    % Initial one is all zeros
        Istr.tstep(t).Mh  = sparse(nsxnjt, nvars);    % Initial one is all zeros
        Istr.tstep(t).Ng  = B2 * G;
        Istr.tstep(t).Nh  = B2 * H;
      else
        % Form D(t)
        D = sparse(nsxnjt, ns*Istr.idx.nj(t-1));
        p1 = Istr.CostWeights(1,1:Istr.idx.nj(t-1),t-1)';
        p1 = p1 / sum(p1);      % sigma(t)
        p2 = Istr.tstep(t).TransMat * p1;
        Di = spdiags(1./p2, 0, Istr.idx.nj(t), Istr.idx.nj(t)) * ...
                sparse(Istr.tstep(t).TransMat) * ...
                spdiags(p1, 0, Istr.idx.nj(t-1), Istr.idx.nj(t-1));
        for i = 1:ns
          D((i-1)*Istr.idx.nj(t)+1:i*Istr.idx.nj(t), (i-1)*Istr.idx.nj(t-1)+1:i*Istr.idx.nj(t-1)) = Di;
        end
        % Apply recursion, form Li, Lf, Mg, Mh, Ng, Nh
        Istr.tstep(t).Li = D  * Istr.tstep(t-1).Lf;
        Istr.tstep(t).Lf = B1 * Istr.tstep(t).Li;
        Istr.tstep(t).Mg = D * Istr.tstep(t-1).Ng;
        Istr.tstep(t).Mh = D * Istr.tstep(t-1).Nh;
        Istr.tstep(t).Ng = B1 * Istr.tstep(t).Mg + B2 * G;
        Istr.tstep(t).Nh = B1 * Istr.tstep(t).Mh + B2 * H;
      end
      Istr.tstep(t).G = G;
      Istr.tstep(t).H = H;
    end
  end

  % Now for the constraint indexing and creation.
  if verbose
    fprintf('- Building constraint submatrices.\n');
  end
  baseMVA = Istr.mpc.baseMVA;
  om = add_constraints(om, 'Pmis', {nt, nj_max, nc_max+1});
  if Istr.DCMODEL
    % Construct all load flow equations using a DC flow model
    if verbose
      fprintf('  - Building DC flow constraints.\n');
    end
    om = add_constraints(om, 'Pf', {nt, nj_max, nc_max+1});
    for t = 1:nt
      for j = 1:Istr.idx.nj(t)
        for k = 1:Istr.idx.nc(t,j)+1
          % First the flow constraints
          mpc = Istr.flow(t,j,k).mpc;
          ion = find(mpc.branch(:, BR_STATUS));
          [Bdc, Bl, Psh, PLsh] = makeBdc(baseMVA, mpc.bus, mpc.branch(ion,:));
          Istr.flow(t,j,k).PLsh = PLsh;     %% save for computing flows later
          negCg = sparse(mpc.gen(:,GEN_BUS), (1:ng)', -1, ...
                        Istr.idx.nb(t,j,k), ng);
          A = [Bdc negCg];
          b = -(mpc.bus(:,PD)+mpc.bus(:,GS))/baseMVA-Psh;
          vs = struct('name', {'Va', 'Pg'}, 'idx', {{t,j,k}, {t,j,k}});
          om = add_constraints(om, 'Pmis', {t,j,k}, A, b, b, vs);
          % Then the thermal limits
          tmp = mpc.branch(ion,RATE_A)/baseMVA;
          iuncon = find(~tmp);
          tmp(iuncon) = Inf(size(iuncon));
          vs = struct('name', {'Va'}, 'idx', {{t,j,k}});
          om = add_constraints(om, 'Pf', {t,j,k}, Bl, -tmp-PLsh, tmp-PLsh, vs);
        end
      end
    end
  else
    if verbose
      fprintf('  - Building load balance constraints.\n');
    end
    % Set simple generation - demand = 0 equations, one for each flow
    for t = 1:nt
      for j = 1:Istr.idx.nj(t)
        for k = 1:Istr.idx.nc(t,j)+1
          mpc = Istr.flow(t,j,k).mpc;
          A = sparse(ones(1, ng));
          b = 1.0*sum(mpc.bus(:, PD)+mpc.bus(:,GS))/baseMVA;
          vs = struct('name', {'Pg'}, 'idx', {{t,j,k}});
          om = add_constraints(om, 'Pmis', {t,j,k}, A, b, b, vs);
        end
      end
    end
  end
  if Istr.IncludeFixedReserves
    if verbose
      fprintf('  - Building fixed zonal reserve constraints.\n');
    end
    om = add_constraints(om, 'Pg_plus_R', {nt, nj_max, nc_max+1});
    om = add_constraints(om, 'Rreq', {nt, nj_max, nc_max+1});
    for t = 1:nt
      for j = 1:Istr.idx.nj(t)
        for k = 1:Istr.idx.nc(t,j)+1
          % First the flow constraints
          mpc = Istr.flow(t,j,k).mpc;
          r = Istr.FixedReserves(t,j,k);
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
          om = add_constraints(om, 'Pg_plus_R', {t,j,k}, A, [], u, vs);
          A = r.zones(:, r.igr);
          l = r.req / mpc.baseMVA;
          vs = struct('name', {'R'}, 'idx', {{t,j,k}});
          om = add_constraints(om, 'Rreq', {t,j,k}, A, l, [], vs);
        end
      end
    end
  end
  
  % Set relationships between generator injections and charge/discharge
  % variables (-pg + psc + psd = 0)
  if verbose && ~isempty(Istr.Storage.UnitIdx)
    fprintf('  - Splitting storage injections into charge/discharge.\n');
  end
  om = add_constraints(om, 'Ps', {nt, nj_max, nc_max+1});
  for t = 1:nt
    for j = 1:Istr.idx.nj(t)
      for k = 1:Istr.idx.nc(t,j)+1
        A = [sparse((1:ns)', Istr.Storage.UnitIdx, -1, ns, ng) Ins Ins];
        b = zeros(ns, 1);
        vs = struct('name', {'Pg', 'Psc', 'Psd'}, 'idx', {{t,j,k}, {t,j,k}, {t,j,k}});
        om = add_constraints(om, 'Ps', {t,j,k}, A, b, b, vs);
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
  om = add_constraints(om, 'ycon', {nt, nj_max, nc_max+1});
  for t = 1:nt,
    for j = 1:Istr.idx.nj(t)
      for k = 1:Istr.idx.nc(t,j)+1
        mpc = Istr.flow(t,j,k).mpc;
        [A, u] = makeAy(baseMVA, ng, mpc.gencost, 1, [], ng+1);
        vs = struct('name', {'Pg', 'y'}, 'idx', {{t,j,k}, {t,j,k}});
        om = add_constraints(om, 'ycon', {t,j,k}, A, [], u, vs);
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
  om = add_constraints(om, 'rampcont', {nt, nj_max, nc_max+1});
  for t =1:nt
    for j = 1:Istr.idx.nj(t)
      for k = 2:Istr.idx.nc(t,j)+1
        ii = find(Istr.flow(t,j,k).mpc.gen(:, GEN_STATUS) > 0);
        ngtmp = length(ii);
        A = sparse((1:ngtmp)', ii, 1, ngtmp, ng);
        u = Istr.flow(t,j,k).mpc.gen(ii,RAMP_10)/baseMVA;
        vs = struct('name', {'Pg', 'Pg'}, 'idx', {{t,j,1}, {t,j,k}});
        om = add_constraints(om, 'rampcont', {t,j,k}, [-A A], -u, u, vs);
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
  om = add_constraints(om, 'dPpRp', {nt, nj_max, nc_max+1});
  for t = 1:nt
    for j = 1:Istr.idx.nj(t);
      for k = 1:Istr.idx.nc(t,j)+1
        ii = find(Istr.flow(t,j,k).mpc.gen(:, GEN_STATUS) > 0);
        ngtmp = length(ii);
        A = sparse((1:ngtmp)', ii, 1, ngtmp, ng);
        l = zeros(ngtmp, 1);
        vs = struct('name', {'dPp', 'Rpp'}, 'idx', {{t,j,k}, {t}});
        om = add_constraints(om, 'dPpRp', {t,j,k}, [-A A], l, [], vs);
      end
    end
  end
  % Negative reserve variables are larger than all decrement variables in
  % all scenarios and flows of a given time slice  0 <= rpm - dpm; these
  % are the ones that set the price of reserves. Include all units that are
  % potentially committed.
  om = add_constraints(om, 'dPmRm', {nt, nj_max, nc_max+1});
  for t = 1:nt
    for j = 1:Istr.idx.nj(t);
      for k = 1:Istr.idx.nc(t,j)+1
        ii = find(Istr.flow(t,j,k).mpc.gen(:, GEN_STATUS) > 0);
        ngtmp = length(ii);
        A = sparse((1:ngtmp)', ii, 1, ngtmp, ng);
        l = zeros(ngtmp, 1);
        vs = struct('name', {'dPm', 'Rpm'}, 'idx', {{t,j,k}, {t}});
        om = add_constraints(om, 'dPmRm', {t,j,k}, [-A A], l, [], vs);
      end
    end
  end
  % The difference between the injection and the contract
  % is equal to the inc minus the dec: Ptjk - Ptc = dPp - dPm
  % Include all units that are potentially committed.
  om = add_constraints(om, 'dPdef', {nt, nj_max, nc_max+1});
  for t = 1:nt
    for j = 1:Istr.idx.nj(t);
      for k = 1:Istr.idx.nc(t,j)+1
        ii = find(Istr.flow(t,j,k).mpc.gen(:, GEN_STATUS) > 0);
        ngtmp = length(ii);
        A = sparse((1:ngtmp)', ii, 1, ngtmp, ng);
        b = zeros(ngtmp, 1);
        vs = struct('name', {'Pg', 'Pc', 'dPp', 'dPm'}, ...
                    'idx', {{t,j,k}, {t}, {t,j,k}, {t,j,k}});
        om = add_constraints(om, 'dPdef', {t,j,k}, [A -A -A A], b, b, vs);
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
  om = add_constraints(om, 'Rrp', {nt, nj_max, nj_max});
  % First, do from t=1:nt-1, since the last one is different and may not
  % even exist depending on the type of horizon
  for t = 1:nt-1
    for j1 = 1:Istr.idx.nj(t)      % j1 is at time t
      for j2 = 1:Istr.idx.nj(t+1)  % j2 is at time t+1
        if Istr.tstep(t+1).TransMask(j2,j1)
          A = [Ing -Ing Ing];
          l = zeros(ng, 1);
          vs = struct('name', {'Pg', 'Pg', 'Rrp'}, ...
                      'idx', {{t,j1,1}, {t+1,j2,1}, {t}});
          om = add_constraints(om, 'Rrp', {t,j1,j2}, A, l, [], vs);
        end
      end
    end
  end
  % Now, pay special attention to a possible last type of ramping
  % constraint. If the horizon involves a terminal value at t=nt+1, then
  % this must also be enforced; in this case, additional ramping
  % reserves procured for t=nt must be defined.  If this
  % condition does not apply, then these reserves are not needed.
  if ~Istr.OpenEnded
    % pterminal <= rrp(nt) + p(nt,j1,0)
    for j1 = 1:Istr.idx.nj(nt)
      A = [Ing Ing];
      l = Istr.TerminalPg/baseMVA;
      vs = struct('name', {'Pg', 'Rrp'}, ...
                  'idx', {{nt,j1,1}, {nt}});
      om = add_constraints(om, 'Rrp', {nt,j1,1}, A, l, [], vs);
    end
  end
  % Now on to downward ramping reserves.
  % First, bound downward ramping reserves from below by all base-case
  % ramping possibilities, 0 <= rrm(t) + p(t+1)j20 - p(t)j10
  om = add_constraints(om, 'Rrm', {nt, nj_max, nj_max});
  % First, do from t=1:nt-1, since the last one is different and may not
  % even exist depending on the type of horizon
  for t = 1:nt-1
    for j1 = 1:Istr.idx.nj(t)      % j1 is at time t
      for j2 = 1:Istr.idx.nj(t+1)  % j2 is at time t+1
        if Istr.tstep(t+1).TransMask(j2,j1)
          A = [-Ing Ing Ing];
          l = zeros(ng, 1);
          vs = struct('name', {'Pg', 'Pg', 'Rrm'}, ...
                      'idx', {{t,j1,1}, {t+1,j2,1}, {t}});
          om = add_constraints(om, 'Rrm', {t,j1,j2}, A, l, [], vs);
        end
      end
    end
  end
  % Now, pay special attention to a possible last type of ramping
  % constraint. If the horizon involves a terminal value at t=nt+1, then
  % this must also be enforced; in this case, additional ramping
  % reserves procured for t=nt must be defined.  If this
  % condition does not apply, then these reserves are not needed.
  if ~Istr.OpenEnded
    % -pterminal <= rrm(nt) - p(nt,j1,0)
    for j1 = 1:Istr.idx.nj(nt)
      A = [-Ing Ing];
      l = -Istr.TerminalPg/baseMVA;
      vs = struct('name', {'Pg', 'Rrm'}, ...
                  'idx', {{nt,j1,1}, {nt}});
      om = add_constraints(om, 'Rrm', {nt,j1,1}, A, l, [], vs);
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
    om = add_constraints(om, 'Sm', {nt, nj_max});
    if Istr.Storage.ForceCyclicStorage
      % sm(1) - beta1*s0 + beta2*Delta_T*[eta_c*psc(1,j,0) + (1/eta_d)*psd(1,j,0)] <= 0
      for j = 1:Istr.idx.nj(1)
        A = [ diagBeta2EtaIn1 diagBeta2overEtaOut1 Ins -spdiags(beta1(:,1), 0, ns, ns)];
        u = zeros(ns, 1);
        vs = struct('name', {'Psc', 'Psd', 'Sm', 'S0'}, 'idx', {{1,j,1}, {1,j,1}, {1}, {}});
        om = add_constraints(om, 'Sm', {1,j}, A, [], u, vs);
      end
    else
      % sm(1) + beta2*Delta_T*[eta_c*psc(1,j,0) + (1/eta_d)*psd(1,j,0)] <= beta1*Initial/baseMVA
      for j = 1:Istr.idx.nj(1)
        A = [ diagBeta2EtaIn1 diagBeta2overEtaOut1 Ins ];
        u = beta1(:,1).*Istr.Storage.InitialStorageLowerBound/baseMVA;
        vs = struct('name', {'Psc', 'Psd', 'Sm'}, 'idx', {{1,j,1}, {1,j,1}, {1}});
        om = add_constraints(om, 'Sm', {1,j}, A, [], u, vs);
      end
    end
    % Then the rest of the periods
    % sm(t) - beta1*(rho(t)*sm(t-1) + (1-rho(t))*s_I(t,j)) + beta2*Delta_T*[eta_c*psc(t,j,0) + (1/eta_d)*psd(t,j,0)] <= 0
    % where s_I(t,j) = L_I(t,j) * s0 + (Mg(t,j)+Mh(t,j)) * x
    for t = 2:nt
      for j = 1:Istr.idx.nj(t)
        Mj = Istr.tstep(t).Mg( j:Istr.idx.nj(t):(ns-1)*Istr.idx.nj(t)+j, :) + ...
             Istr.tstep(t).Mh( j:Istr.idx.nj(t):(ns-1)*Istr.idx.nj(t)+j, :);
        Lij = Istr.tstep(t).Li( j:Istr.idx.nj(t):(ns-1)*Istr.idx.nj(t)+j, :);
        diag1minusRhoBeta1 = spdiags((1-rho(:,t)) .* beta1(:,t), 0, ns, ns);
        A = sparse([1:ns,1:ns,1:ns,1:ns]', ...
                   [vv.i1.Psc(t,j,1):vv.iN.Psc(t,j,1), vv.i1.Psd(t,j,1):vv.iN.Psd(t,j,1), vv.i1.Sm(t-1):vv.iN.Sm(t-1), vv.i1.Sm(t):vv.iN.Sm(t)]', ...
                   [beta2EtaIn(:,t); beta2overEtaOut(:,t); -beta1(:,t).*rho(:,t); ones(ns,1)], ...
                   ns, nvars) ...
                - diag1minusRhoBeta1 * Mj;
        if Istr.Storage.ForceCyclicStorage
          As0 = sparse(ns, nvars);
          As0(:, vv.i1.S0:vv.iN.S0) = -diag1minusRhoBeta1 * Lij;
          A = A + As0;
          u = zeros(ns, 1);
        else
          u = full(diag1minusRhoBeta1 * Lij * Istr.Storage.InitialStorage/baseMVA);
        end
        om = add_constraints(om, 'Sm', {t,j}, A, [], u);
      end
    end
    % Do the same we did for sm(t) for sp(t). First the initial step ...
    om = add_constraints(om, 'Sp', {nt, nj_max});
    if Istr.Storage.ForceCyclicStorage
      % -sp(1) + beta1*s0 - beta2*Delta_T*[eta_c*psc(1,j,0) + (1/eta_d)*psd(1,j,0)] <= 0
      for j = 1:Istr.idx.nj(1)
        A = [ -diagBeta2EtaIn1 -diagBeta2overEtaOut1 -Ins spdiags(beta1(:,1), 0, ns, ns) ];
        u = zeros(ns, 1);
        vs = struct('name', {'Psc', 'Psd', 'Sp', 'S0'}, 'idx', {{1,j,1}, {1,j,1}, {1}, {}});
        om = add_constraints(om, 'Sp', {1,j}, A, [], u, vs);
      end
    else
      % -sp(1) - beta2*Delta_T*[eta_c*psc(1,j,0) + (1/eta_d)*psd(1,j,0)] <= -beta1*Initial/baseMVA
      for j = 1:Istr.idx.nj(1)
        A = [ -diagBeta2EtaIn1 -diagBeta2overEtaOut1 -Ins ];
        u = -beta1(:,1).*Istr.Storage.InitialStorageUpperBound/baseMVA;
        vs = struct('name', {'Psc', 'Psd', 'Sp'}, 'idx', {{1,j,1}, {1,j,1}, {1}});
        om = add_constraints(om, 'Sp', {1,j}, A, [], u, vs);
      end
    end
    % Then the rest of the periods
    % -sp(t) + beta1*(rho(t)*sp(t-1) + (1-rho(t))*s_I(t,j)) - beta2*Delta_T*[eta_c*psc(t,j,0) + (1/eta_d)*psd(t,j,0)] <= 0
    % where s_I(t,j) = L_I(t,j) * s0 + (Mg(t,j)+Mh(t,j)) * x
    for t = 2:nt
      for j = 1:Istr.idx.nj(t)
        Mj = Istr.tstep(t).Mg( j:Istr.idx.nj(t):(ns-1)*Istr.idx.nj(t)+j, :) + ...
             Istr.tstep(t).Mh( j:Istr.idx.nj(t):(ns-1)*Istr.idx.nj(t)+j, :);
        Lij = Istr.tstep(t).Li( j:Istr.idx.nj(t):(ns-1)*Istr.idx.nj(t)+j, :);
        diag1minusRhoBeta1 = spdiags((1-rho(:,t)) .* beta1(:,t), 0, ns, ns);
        A = sparse([1:ns,1:ns,1:ns,1:ns]', ...
                   [vv.i1.Psc(t,j,1):vv.iN.Psc(t,j,1), vv.i1.Psd(t,j,1):vv.iN.Psd(t,j,1), vv.i1.Sp(t-1):vv.iN.Sp(t-1), vv.i1.Sp(t):vv.iN.Sp(t)]', ...
                   [-beta2EtaIn(:,t); -beta2overEtaOut(:,t); beta1(:,t).*rho(:,t); -ones(ns,1)], ...
                   ns, nvars) ...
                + diag1minusRhoBeta1 * Mj;
        if Istr.Storage.ForceCyclicStorage
          As0 = sparse(ns, nvars);
          As0(:, vv.i1.S0:vv.iN.S0) = diag1minusRhoBeta1 * Lij;
          A = A + As0;
          u = zeros(ns, 1);
        else
          u = full(-diag1minusRhoBeta1 * Lij * Istr.Storage.InitialStorage/baseMVA);
        end
        om = add_constraints(om, 'Sp', {t,j}, A, [], u);
      end
    end
    % Now go on and limit the amount of energy that can be used if a
    % contingency does happen. Bound sm first. First examine time period 1 wrt to initial
    % stored energy, then t=2 and on.
    om = add_constraints(om, 'contSm', {nt, nj_max, nc_max+1});
    for j = 1:Istr.idx.nj(1)
      for k = 2:Istr.idx.nc(1,j)+1  %% NOTE NO k=1!!!
        if Istr.Storage.ForceCyclicStorage
          % beta4*Delta_T*[eta_c*psc(1,j,0) + (1/eta_d)*psd(1,j,0)] + beta3*Delta_T*[eta_c*psc(1,j,k) + (1/eta_d)*psd(1,j,k)] - beta5*s0 <= -sm_min(1)
          A = [ diagBeta4EtaIn1 diagBeta4overEtaOut1 diagBeta3EtaIn1 diagBeta3overEtaOut1 -spdiags(beta5(:,1), 0, ns, ns) ];
          u = -MinStorageLevel(:,1)/baseMVA;
          vs = struct('name', {'Psc', 'Psd', 'Psc', 'Psd', 'S0'}, ...
                      'idx', {{1,j,1}, {1,j,1}, {1,j,k}, {1,j,k}, {}});
        else
          % beta4*Delta_T*[eta_c*psc(1,j,0) + (1/eta_d)*psd(1,j,0)] + beta3*Delta_T*[eta_c*psc(1,j,k) + (1/eta_d)*psd(1,j,k)] <= beta5*Initial/baseMVA - sm_min(1)
          A = [ diagBeta4EtaIn1 diagBeta4overEtaOut1 diagBeta3EtaIn1 diagBeta3overEtaOut1 ];
          u = (beta5(:,1).*Istr.Storage.InitialStorageLowerBound - MinStorageLevel(:,1))/baseMVA;
          vs = struct('name', {'Psc', 'Psd', 'Psc', 'Psd'}, ...
                      'idx', {{1,j,1}, {1,j,1}, {1,j,k}, {1,j,k}});
        end
        om = add_constraints(om, 'contSm', {1,j,k}, A, [], u, vs);
      end
    end
    % then the rest of the periods
    % -beta5*(rho(t)*sm(t-1) + (1-rho(t))*s_I(t,j)) + beta4*Delta_T*[eta_c*psc(t,j,0) + (1/eta_d)*psd(t,j,0)] + beta3*Delta_T*[eta_c*psc(t,j,k) + (1/eta_d)*psd(t,j,k)] <= -sm_min(t)
    % where s_I(t,j) = L_I(t,j) * s0 + (Mg(t,j)+Mh(t,j)) * x
    for t = 2:nt
      for j = 1:Istr.idx.nj(t)
        for k = 2:Istr.idx.nc(t,j)+1
          Mj = Istr.tstep(t).Mg( j:Istr.idx.nj(t):(ns-1)*Istr.idx.nj(t)+j, :) + ...
               Istr.tstep(t).Mh( j:Istr.idx.nj(t):(ns-1)*Istr.idx.nj(t)+j, :);
          Lij = Istr.tstep(t).Li( j:Istr.idx.nj(t):(ns-1)*Istr.idx.nj(t)+j, :);
          diag1minusRhoBeta5 = spdiags((1-rho(:,t)) .* beta5(:,t), 0, ns, ns);
          A = sparse([1:ns,1:ns,1:ns,1:ns,1:ns]', ...
                     [vv.i1.Psc(t,j,1):vv.iN.Psc(t,j,1), vv.i1.Psd(t,j,1):vv.iN.Psd(t,j,1), vv.i1.Psc(t,j,k):vv.iN.Psc(t,j,k), vv.i1.Psd(t,j,k):vv.iN.Psd(t,j,k), vv.i1.Sm(t-1):vv.iN.Sm(t-1)]', ...
                     [beta4EtaIn(:,t); beta4overEtaOut(:,t); beta3EtaIn(:,t); beta3overEtaOut(:,t); -beta5(:,t).*rho(:,t)], ...
                     ns, nvars) ...
                  - diag1minusRhoBeta5 * Mj;
          u = -MinStorageLevel(:,t)/baseMVA;
          if Istr.Storage.ForceCyclicStorage
            As0 = sparse(ns, nvars);
            As0(:, vv.i1.S0:vv.iN.S0) = -diag1minusRhoBeta5 * Lij;
            A = A + As0;
          else
            u = u + diag1minusRhoBeta5 * Lij * Istr.Storage.InitialStorageLowerBound/baseMVA;
          end
          om = add_constraints(om, 'contSm', {t,j,k}, A, [], u);
        end
      end
    end
    % Bound sp first. First examine time period 1 wrt to initial
    % stored energy, then t=2 and on.
    om = add_constraints(om, 'contSp', {nt, nj_max, nc_max+1});
    for j = 1:Istr.idx.nj(1)
      for k = 2:Istr.idx.nc(1,j)+1
        if Istr.Storage.ForceCyclicStorage
          % -beta4*Delta_T*[eta_c*psc(1,j,0) + (1/eta_d)*psd(1,j,0)] - beta3*Delta_T*[eta_c*psc(1,j,k) + (1/eta_d)*psd(1,j,k)] + beta5*s0 <= sp_max(1)
          A = [ -diagBeta4EtaIn1 -diagBeta4overEtaOut1 -diagBeta3EtaIn1 -diagBeta3overEtaOut1 spdiags(beta5(:,1), 0, ns, ns)];
          u = MaxStorageLevel(:,1)/baseMVA;
          vs = struct('name', {'Psc', 'Psd', 'Psc', 'Psd', 'S0'}, ...
                      'idx', {{1,j,1}, {1,j,1}, {1,j,k}, {1,j,k}, {}});
        else
          % -beta4*Delta_T*[eta_c*psc(1,j,0) + (1/eta_d)*psd(1,j,0)] - beta3*Delta_T*[eta_c*psc(1,j,k) + (1/eta_d)*psd(1,j,k)] <= -beta5*Initial/baseMVA + sp_max(1)
          A = [ -diagBeta4EtaIn1 -diagBeta4overEtaOut1 -diagBeta3EtaIn1 -diagBeta3overEtaOut1 ];
          u = (MaxStorageLevel(:,1) - beta5(:,1).*Istr.Storage.InitialStorageUpperBound)/baseMVA;
          vs = struct('name', {'Psc', 'Psd', 'Psc', 'Psd'}, ...
                      'idx', {{1,j,1}, {1,j,1}, {1,j,k}, {1,j,k}});
        end
        om = add_constraints(om, 'contSp', {1,j,k}, A, [], u, vs);
      end
    end
    % then the rest of the periods
    % beta5*(rho(t)*sp(t-1) + (1-rho(t))*s_I(t,j)) - beta4*Delta_T*[eta_c*psc(t,j,0) + (1/eta_d)*psd(t,j,0)] - beta3*Delta_T*[eta_c*psc(t,j,k) + (1/eta_d)*psd(t,j,k)] <= sp_max(t)
    % where s_I(t,j) = L_I(t,j) * s0 + (Mg(t,j)+Mh(t,j)) * x
    for t = 2:nt
      for j = 1:Istr.idx.nj(t)
        for k = 2:Istr.idx.nc(t,j)+1
          Mj = Istr.tstep(t).Mg( j:Istr.idx.nj(t):(ns-1)*Istr.idx.nj(t)+j, :) + ...
               Istr.tstep(t).Mh( j:Istr.idx.nj(t):(ns-1)*Istr.idx.nj(t)+j, :);
          Lij = Istr.tstep(t).Li( j:Istr.idx.nj(t):(ns-1)*Istr.idx.nj(t)+j, :);
          diag1minusRhoBeta5 = spdiags((1-rho(:,t)) .* beta5(:,t), 0, ns, ns);
          A = sparse([1:ns,1:ns,1:ns,1:ns,1:ns]', ...
                     [vv.i1.Psc(t,j,1):vv.iN.Psc(t,j,1), vv.i1.Psd(t,j,1):vv.iN.Psd(t,j,1), vv.i1.Psc(t,j,k):vv.iN.Psc(t,j,k), vv.i1.Psd(t,j,k):vv.iN.Psd(t,j,k), vv.i1.Sp(t-1):vv.iN.Sp(t-1)]', ...
                     [-beta4EtaIn(:,t); -beta4overEtaOut(:,t); -beta3EtaIn(:,t); -beta3overEtaOut(:,t); beta5(:,t).*rho(:,t)], ...
                     ns, nvars) ...
                  + diag1minusRhoBeta5 * Mj;
          u = MaxStorageLevel(:,t)/baseMVA;
          if Istr.Storage.ForceCyclicStorage
            As0 = sparse(ns, nvars);
            As0(:, vv.i1.S0:vv.iN.S0) = diag1minusRhoBeta5 * Lij;
            A = A + As0;
          else
            u = u - diag1minusRhoBeta5 * Lij * Istr.Storage.InitialStorageUpperBound/baseMVA;
          end
          om = add_constraints(om, 'contSp', {t,j,k}, A, [], u);
        end
      end
    end
  end
  
  % Now, if required, constrain the expected terminal storage quantity; two
  % different ways:
  if Istr.Storage.ForceExpectedTerminalStorage && Istr.Storage.ForceCyclicStorage
    error('most: ForceExpectedTerminalStorage and ForceCyclicStorage cannot be simultaneously true.');
  end
  if ns
    % The following code assumes that no more variables will be added
    if Istr.Storage.ForceExpectedTerminalStorage
      % 1) Constrain the expected terminal storage to be some target value
      A = sparse(ns, nvars);
      b = zeros(ns, 1);
      for j = 1:Istr.idx.nj(nt)
        Ngj = Istr.tstep(nt).Ng( j:Istr.idx.nj(nt):(ns-1)*Istr.idx.nj(nt)+j, :);
        Nhj = Istr.tstep(nt).Nh( j:Istr.idx.nj(nt):(ns-1)*Istr.idx.nj(nt)+j, :);
        Lfj = Istr.tstep(nt).Lf( j:Istr.idx.nj(nt):(ns-1)*Istr.idx.nj(nt)+j, :);
        A = A + Istr.CostWeights(1,j,nt) * (Ngj + Nhj);
        b = b + Istr.CostWeights(1,j,nt) * (Lfj * Istr.Storage.InitialStorage) / baseMVA;
      end
      endprob = sum(Istr.CostWeights(1,1:Istr.idx.nj(nt),nt)');
      A = (1/endprob) * A;
      b = (1/endprob) * b;
      l = Istr.Storage.ExpectedTerminalStorageMin / baseMVA - b;
      u = Istr.Storage.ExpectedTerminalStorageMax / baseMVA - b;
      om = add_constraints(om, 'ESnt', A, l, u);
    elseif Istr.Storage.ForceCyclicStorage
      % 2) Constrain the initial storage (a variable) to be the same as the final expected storage
      A = sparse(ns, nvars);
      for j = 1:Istr.idx.nj(nt)
        Ngj = Istr.tstep(nt).Ng( j:Istr.idx.nj(nt):(ns-1)*Istr.idx.nj(nt)+j, :);
        Nhj = Istr.tstep(nt).Nh( j:Istr.idx.nj(nt):(ns-1)*Istr.idx.nj(nt)+j, :);
        A = A + Istr.CostWeights(1,j,nt) * (Ngj + Nhj);
      end
      endprob = sum(Istr.CostWeights(1,1:Istr.idx.nj(nt),nt)');
      A = (1/endprob) * A;
      for j = 1:Istr.idx.nj(nt)
        Lfj = Istr.tstep(nt).Lf(j:Istr.idx.nj(nt):(ns-1)*Istr.idx.nj(nt)+j, :);
        A(:, vv.i1.S0:vv.iN.S0) = A(:, vv.i1.S0:vv.iN.S0) ...
                  + Istr.CostWeights(1,j,nt) * Lfj;
      end
      A(:, vv.i1.S0:vv.iN.S0) = (1/endprob) * A(:, vv.i1.S0:vv.iN.S0) - speye(ns);
      b = zeros(ns, 1);
      om = add_constraints(om, 'ESnt', A, b, b);
    end
  end

  % Dynamical system contraints
  if nzd || nyo
    if verbose
      fprintf('  - Building dynamic system constraints.\n');
    end
    % Compute matrices that give the expected dispatch in time period t
    % given that we make it to that period, for all generators at once,
    % when multiplied by x, i.e. E[p(t)] = E(t) * x
    for t = 1:nt
      E = sparse(ng, nvars);
      for j = 1:Istr.idx.nj(t)
        for k = 1:Istr.idx.nc(t,j)+1;
          E = E + (Istr.CostWeightsAdj(k,j,t)/Istr.StepProb(t)) * sparse((1:ng)', (vv.i1.Pg(t,j,k):vv.iN.Pg(t,j,k))', 1, ng, nvars);
        end
      end
      Istr.tstep(t).E = E;
    end
  end

  % Form the dynamical system state equations and bound constraints on the
  % state vector
  if nzd
    om = add_constraints(om, 'DSz', {nyt-1});
    b = zeros(nzd, 1);
    for t = 1:nyt-1
      if t <= nt  % We have p(t) available to drive the dynamical system up to t=nt
        % Form the constraint matrix, so B*E*x + A*z(t) - I*z(t+1) = 0
        A = Istr.dstep(t).B * Istr.tstep(t).E;
      else
        % The dynamical system horizon is longer than the injection planning
        % horizon and we don't know what p(t) is, but continue to drive the
        % dynamical system as if p(t) = 0 and perhaps take that into account
        % when setting Ymax, Ymin in this time window.  That is, A*z(t) - I*z(t+1) = 0
        A = sparse(nzd, nvars);
      end
      A(:, vv.i1.Z(t):vv.iN.Z(t)) = Istr.dstep(t).A;
      A(:, vv.i1.Z(t+1):vv.iN.Z(t+1)) = -speye(nzd);
      om = add_constraints(om, 'DSz', {t}, A, b, b);
    end
  end
  
  % Form the output equations and their restrictions
  if nyo
    om = add_constraints(om, 'DSy', {nyt});
    for t = 1:nyt
      if t <= nt
        A = Istr.dstep(t).D * Istr.tstep(t).E;
      else
        A = sparse(nyo, nvars);
      end
      if nzd
        A(:, vv.i1.Z(t):vv.iN.Z(t)) = Istr.dstep(t).C;
      end
      l = Istr.dstep(t).ymin;
      u = Istr.dstep(t).ymax;
      om = add_constraints(om, 'DSy', {t}, A, l, u);
    end
  end
  
  % UNIT COMMITMENT
  if UC
    if verbose
      fprintf('  - Building unit commitment constraints.\n');
    end
    % u(t,i) - u(t-1,i) - v(t,i) + w(t,i) = 0
    om = add_constraints(om, 'uvw', {nt});
    for t = 1:nt
      if t == 1
        % First for t=1 when u(t-1,i) is really u(0,i) or u(nt,i)
        if Istr.UC.CyclicCommitment
          vs = struct('name', {'u', 'u', 'v', 'w'}, 'idx', {{1}, {nt}, {1}, {1}});
          A = [Ing -Ing -Ing Ing];
          b = zeros(ng, 1);
        else
          vs = struct('name', {'u', 'v', 'w'}, 'idx', {{1}, {1}, {1}});
          A = [Ing -Ing Ing];
          b = (Istr.UC.InitialState > 0);
        end
      else
        % Then for rest of periods
        vs = struct('name', {'u', 'u', 'v', 'w'}, 'idx', {{t}, {t-1}, {t}, {t}});
        A = [Ing -Ing -Ing Ing];
        b = zeros(ng, 1);
      end
      om = add_constraints(om, 'uvw', {t}, A, b, b, vs);
    end
    % Then continue with minimimum up time constraints. Again, two
    % different forms depending on whether the horizon is cyclical or not
    om = add_constraints(om, 'minup', {nt, ng});
    for t = 1:nt
      for i = 1:ng
        ti = t-Istr.UC.MinUp(i)+1:t;
        if Istr.UC.CyclicCommitment     % window is circular
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
        om = add_constraints(om, 'minup', {t, i}, A, [], 0, vs);
      end
    end
    % Continue with minimimum downtime constraints. Two
    % different forms depending on whether the horizon is cyclical or not
    om = add_constraints(om, 'mindown', {nt, ng});
    for t = 1:nt
      for i = 1:ng
        ti = t-Istr.UC.MinDown(i)+1:t;
        if Istr.UC.CyclicCommitment     % window is circular
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
        om = add_constraints(om, 'mindown', {t, i}, A, [], 1, vs);
      end
    end
    % Limit generation ranges based on commitment status; first Pmax;
    % p - u*Pmax <= 0
    % For contingent flows, however, if a generator is ousted as a result
    % of the contingency, then this constraint should not be enforced.
    om = add_constraints(om, 'uPmax', {nt, nj_max, nc_max+1});
    for t = 1:nt
      for j = 1:Istr.idx.nj(t)
        for k = 1:Istr.idx.nc(t,j)+1
          mpc = Istr.flow(t,j,k).mpc;
          ii = find(mpc.gen(:, GEN_STATUS));
          nii = length(ii);
          vs = struct('name', {'Pg', 'u'}, 'idx', {{t,j,k}, {t}});
          A = [ sparse(1:nii, ii, 1, nii, ng) ...
                sparse(1:nii, ii, -mpc.gen(ii, PMAX)/baseMVA, nii, ng) ];
          u = zeros(nii, 1);
          om = add_constraints(om, 'uPmax', {t,j,k}, A, [], u, vs);
        end
      end
    end
    % Then Pmin,  -p + u*Pmin <= 0
    om = add_constraints(om, 'uPmin', {nt, nj_max, nc_max+1});
    for t = 1:nt
      for j = 1:Istr.idx.nj(t)
        for k = 1:Istr.idx.nc(t,j)+1
          mpc = Istr.flow(t,j,k).mpc;
          ii = find(mpc.gen(:, GEN_STATUS));
          nii = length(ii);
          vs = struct('name', {'Pg', 'u'}, 'idx', {{t,j,k}, {t}});
          A = [ sparse(1:nii, ii, -1, nii, ng) ...
                sparse(1:nii, ii, mpc.gen(ii, PMIN)/baseMVA, nii, ng) ];
          u = zeros(nii, 1);
          om = add_constraints(om, 'uPmin', {t,j,k}, A, [], u, vs);
        end
      end
    end
    % Then, if there is Qg coordination, do the same for Qg
    % q - u*Qmax <= 0
    % For contingent flows, however, if a generator is ousted as a result
    % of the contingency, then this constraint should not be enforced.
    if Istr.QCoordination
      om = add_constraints(om, 'uQmax', {nt, nj_max, nc_max+1});
      for t = 1:nt
        for j = 1:Istr.idx.nj(t)
          for k = 1:Istr.idx.nc(t,j)+1
            mpc = Istr.flow(t,j,k).mpc;
            ii = find(mpc.gen(:, GEN_STATUS));
            nii = length(ii);
            vs = struct('name', {'Qg', 'u'}, 'idx', {{t,j,k}, {t}});
            A = [ sparse(1:nii, ii, 1, nii, ng) ...
                  sparse(1:nii, ii, -mpc.gen(ii, QMAX)/baseMVA, nii, ng) ];
            u = zeros(nii, 1);
            om = add_constraints(om, 'uQmax', {t,j,k}, A, [], u, vs);
          end
        end
      end
      % Then Qmin,  -q + u*Qmin <= 0
      om = add_constraints(om, 'uQmin', {nt, nj_max, nc_max+1});
      for t = 1:nt
        for j = 1:Istr.idx.nj(t)
          for k = 1:Istr.idx.nc(t,j)+1
            mpc = Istr.flow(t,j,k).mpc;
            ii = find(mpc.gen(:, GEN_STATUS));
            nii = length(ii);
            vs = struct('name', {'Qg', 'u'}, 'idx', {{t,j,k}, {t}});
            A = [ sparse(1:nii, ii, -1, nii, ng) ...
                  sparse(1:nii, ii, mpc.gen(ii, QMIN)/baseMVA, nii, ng) ];
            u = zeros(nii, 1);
            om = add_constraints(om, 'uQmin', {t,j,k}, A, [], u, vs);
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
  c1 = 0;

  % First assign the ramping costs; H1 has few coefficients initially and
  % this should make the shuffling and reordering of coefficients more
  % efficient.  All other accesses to H1 will be diagonal insertions, which
  % take less time than anti-diagonal insertions.
  % First do first period wrt to InitialPg.
  om = add_costs(om, 'RampWear', {nt+1, nj_max, nj_max});
  for j = 1:Istr.idx.nj(1)
    w = Istr.tstep(1).TransMat(j,1);  % the probability of going from initial state to jth
    H = spdiags(w * baseMVA^2 * Istr.RampWearCostCoeff(:,1), 0, ng, ng);
    Cw = -w * baseMVA * Istr.RampWearCostCoeff(:,1) .* Istr.InitialPg;
    cp = struct('H', H, 'Cw', Cw);
    vs = struct('name', {'Pg'}, 'idx', {{1,j,1}});
    om = add_costs(om, 'RampWear', {1,j,1}, cp, vs);
    c1 = c1 + w * 0.5 * Istr.RampWearCostCoeff(:,1)' * Istr.InitialPg.^2;
  end
  % Then the remaining periods
  for t = 2:nt
    for j2 = 1:Istr.idx.nj(t)
      for j1 = 1:Istr.idx.nj(t-1)
        w = Istr.tstep(t).TransMat(j2,j1) * Istr.CostWeights(1, j1, t-1);
        h = w * baseMVA^2 * Istr.RampWearCostCoeff(:,t);
        i = (1:ng)';
        j = ng+(1:ng)';
        H = sparse([i;j;i;j], [i;i;j;j], [h;-h;-h;h], 2*ng, 2*ng);
        cp = struct('H', H, 'Cw', zeros(2*ng,1));
        vs = struct('name', {'Pg', 'Pg'}, 'idx', {{t-1,j1,1}, {t,j2,1}});
        om = add_costs(om, 'RampWear', {t,j1,j2}, cp, vs);
      end
    end
  end
  % Finally, if there is a terminal state problem, apply cost to
  % the transition starting from t=nt.  Note that in this case
  % Istr.tstep(nt+1).TransMat must be defined! it is the only piece of data
  % that makes sense for nt+1; all other fields in Istr.tstep(nt+1) can be empty.
  if ~Istr.OpenEnded
    for j = 1:Istr.idx.nj(nt)
      w = Istr.tstep(nt+1).TransMat(1, j) * Istr.CostWeights(1, j, nt);
      H = spdiags(w * baseMVA^2 * Istr.RampWearCostCoeff(:,nt+1), 0, ng, ng);
      Cw = -w * baseMVA * Istr.RampWearCostCoeff(:,nt+1) .* Istr.TerminalPg;
      cp = struct('H', H, 'Cw', Cw);
      vs = struct('name', {'Pg'}, 'idx', {{nt,j,1}});
      om = add_costs(om, 'RampWear', {nt+1,j,1}, cp, vs);
      c1 = c1 + w * 0.5 * Istr.RampWearCostCoeff(:,nt+1)' * Istr.TerminalPg.^2;
    end
  end

  % Now go on and assign energy, inc/dec and contingency reserves
  % costs for all committed units.
  om = add_costs(om, 'Cp', {nt, nj_max, nc_max+1});
  om = add_costs(om, 'Cy', {nt, nj_max, nc_max+1});
  om = add_costs(om, 'Cpp', {nt, nj_max, nc_max+1});
  om = add_costs(om, 'Cpm', {nt, nj_max, nc_max+1});
  if Istr.IncludeFixedReserves
    om = add_costs(om, 'Rcost', {nt, nj_max, nc_max+1});
  end
  om = add_costs(om, 'Crpp', {nt});
  om = add_costs(om, 'Crpm', {nt});
  for t = 1:nt
    for j = 1:Istr.idx.nj(t)
      for k = 1:Istr.idx.nc(t,j)+1
        w = Istr.CostWeightsAdj(k,j,t);     %% NOTE (k,j,t) order !!!

        % weighted polynomial energy costs for committed units
        gc = Istr.flow(t,j,k).mpc.gencost;
        ipol = find(gc(:, MODEL) == POLYNOMIAL);
        if ~isempty(ipol)
          ncost = gc(ipol(1), NCOST);
          if all(gc(ipol, NCOST) == ncost)    %% uniform order of polynomials
            %% use vectorized code
            if ncost > 3
              error('most: polynomial generator costs of order higher than quadratic not supported');
            elseif ncost == 3
              H = sparse(ipol, ipol, 2 * w * baseMVA^2*gc(ipol, COST), ng, ng);
            else
              H = sparse(ng,ng);
            end
            Cw = zeros(ng, 1);
            if ncost >= 2
              Cw(ipol) = w * baseMVA*gc(ipol, COST+ncost-2);
            end
            c1 = c1 + w * sum(gc(ipol, COST+ncost-1));
          else                                %% non-uniform order of polynomials
            %% use a loop
            H = sparse(ng,ng);
            Cw = zeros(ng, 1);
            for i = ipol'
              ncost = gc(i, NCOST);
              if ncost > 3
                error('most: polynomial generator costs of order higher than quadratic not supported');
              elseif ncost == 3
                H(i,i) = 2 * w * baseMVA^2*gc(i, COST);
              end
              if ncost >= 2
                Cw(i) = w * baseMVA*gc(i, COST+ncost-2);
              end
              c1 = c1 + w * gc(i, COST+ncost-1);
            end
          end
          cp = struct('H', H, 'Cw', Cw);
          vs = struct('name', {'Pg'}, 'idx', {{t,j,k}});
          om = add_costs(om, 'Cp', {t,j,k}, cp, vs);
        end

        % weighted y-variables for piecewise linear energy costs for committed units
        % ipwl = find( (Istr.flow(t,j,k).mpc.gen(:,GEN_STATUS) > 0) & (gc(:,MODEL) == PW_LINEAR));
        if Istr.idx.ny(t,j,k)
          cp = struct('Cw', w * ones(Istr.idx.ny(t,j,k),1));
          vs = struct('name', {'y'}, 'idx', {{t,j,k}});
          om = add_costs(om, 'Cy', {t,j,k}, cp, vs);
        end

        % inc and dec offers for each flow
        cp = struct('Cw', w * baseMVA * Istr.offer(t).PositiveActiveDeltaPrice(:));
        vs = struct('name', {'dPp'}, 'idx', {{t,j,k}});
        om = add_costs(om, 'Cpp', {t,j,k}, cp, vs);
        cp = struct('Cw', w * baseMVA * Istr.offer(t).NegativeActiveDeltaPrice(:));
        vs = struct('name', {'dPm'}, 'idx', {{t,j,k}});
        om = add_costs(om, 'Cpm', {t,j,k}, cp, vs);

        % weighted fixed reserves cost
        if Istr.IncludeFixedReserves
          cp = struct('Cw', w * Istr.FixedReserves(t,j,k).cost(r.igr) * baseMVA);
          vs = struct('name', {'R'}, 'idx', {{t,j,k}});
          om = add_costs(om, 'Rcost', {t,j,k}, cp, vs);
        end
      end
    end
    
    % contingency reserve costs
    cp = struct('Cw', baseMVA * Istr.StepProb(t) * Istr.offer(t).PositiveActiveReservePrice(:));
    vs = struct('name', {'Rpp'}, 'idx', {{t}});
    om = add_costs(om, 'Crpp', {t}, cp, vs);
    cp = struct('Cw', baseMVA * Istr.StepProb(t) * Istr.offer(t).NegativeActiveReservePrice(:));
    vs = struct('name', {'Rpm'}, 'idx', {{t}});
    om = add_costs(om, 'Crpm', {t}, cp, vs);
  end
  % Assign load following ramp reserve costs.  Do first nt-1 periods first
  om = add_costs(om, 'Crrp', {Istr.idx.ntramp});
  om = add_costs(om, 'Crrm', {Istr.idx.ntramp});
  for t = 1:nt-1,
    cp = struct('Cw', baseMVA * Istr.StepProb(t+1) * Istr.offer(t).PositiveLoadFollowReservePrice(:));
    vs = struct('name', {'Rrp'}, 'idx', {{t}});
    om = add_costs(om, 'Crrp', {t}, cp, vs);
    cp = struct('Cw', baseMVA * Istr.StepProb(t+1) * Istr.offer(t).NegativeLoadFollowReservePrice(:));
    vs = struct('name', {'Rrm'}, 'idx', {{t}});
    om = add_costs(om, 'Crrm', {t}, cp, vs);
  end
  % Then do last period if needed Terminal state case
  if ~Istr.OpenEnded
    %% are these costs missing a Istr.StepProb(t)?  -- rdz
    cp = struct('Cw', baseMVA * Istr.offer(nt).PositiveLoadFollowReservePrice(:));
    vs = struct('name', {'Rrp'}, 'idx', {{nt}});
    om = add_costs(om, 'Crrp', {nt}, cp, vs);
    cp = struct('Cw', baseMVA * Istr.offer(nt).NegativeLoadFollowReservePrice(:));
    vs = struct('name', {'Rrm'}, 'idx', {{nt}});
    om = add_costs(om, 'Crrm', {nt}, cp, vs);
  end
  % Assign startup/shutdown costs, if any, and fixed operating costs
  if UC
    om = add_costs(om, 'c00', {nt});
    om = add_costs(om, 'startup', {nt});
    om = add_costs(om, 'shutdown', {nt});
    for t = 1:nt
      ww = zeros(ng, 1);
      for j = 1:Istr.idx.nj(t)
        for k = 1:Istr.idx.nc(t)+1
          ww = ww + Istr.CostWeightsAdj(k,j,t) * Istr.flow(t,j,k).mpc.gen(:, GEN_STATUS);
        end
      end
      cp = struct('Cw', ww.*Istr.UC.c00(:,t));
      vs = struct('name', {'u'}, 'idx', {{t}});
      om = add_costs(om, 'c00', {t}, cp, vs);
      cp = struct('Cw', Istr.StepProb(t)*Istr.flow(t,1,1).mpc.gencost(:, STARTUP));
      vs = struct('name', {'v'}, 'idx', {{t}});
      om = add_costs(om, 'startup', {t}, cp, vs);
      cp = struct('Cw', Istr.StepProb(t)*Istr.flow(t,1,1).mpc.gencost(:, SHUTDOWN));
      vs = struct('name', {'w'}, 'idx', {{t}});
      om = add_costs(om, 'shutdown', {t}, cp, vs);
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
    vv = get_idx(om);
    for t = 1:nt
      % Compute cost coefficients for value of expected leftover storage
      % after a contingency
      for j = 1:Istr.idx.nj(t)
        % pick rows for jth base injections
        Gtj0  = Istr.tstep(t).G(  j:Istr.idx.nj(t):(ns-1)*Istr.idx.nj(t)+j, :);
        Htj0  = Istr.tstep(t).H(  j:Istr.idx.nj(t):(ns-1)*Istr.idx.nj(t)+j, :);
        Litj0 = Istr.tstep(t).Li( j:Istr.idx.nj(t):(ns-1)*Istr.idx.nj(t)+j, :);
        Mgtj0 = Istr.tstep(t).Mg( j:Istr.idx.nj(t):(ns-1)*Istr.idx.nj(t)+j, :);
        Mhtj0 = Istr.tstep(t).Mh( j:Istr.idx.nj(t):(ns-1)*Istr.idx.nj(t)+j, :);
        sum_psi_tjk = sum(Istr.CostWeights(2:Istr.idx.nc(t,j)+1,j,t));
        if t == nt
          A1 = A1 + Istr.CostWeights(1,j,t) * spdiags(OutEff(:,t) .* beta1(:,t), 0, ns, ns) * Litj0;
          A2 = A2 + Istr.CostWeights(1,j,t) * spdiags(OutEff(:,t) .* beta1(:,t), 0, ns, ns) * Mgtj0;
          A3 = A3 + Istr.CostWeights(1,j,t) * spdiags(OutEff(:,t) .* beta1(:,t), 0, ns, ns) * Mhtj0;
          A4 = A4 + Istr.CostWeights(1,j,t) * spdiags(OutEff(:,t) .* beta2(:,t), 0, ns, ns) * Gtj0;
          A5 = A5 + Istr.CostWeights(1,j,t) * spdiags(OutEff(:,t) .* beta2(:,t), 0, ns, ns) * Htj0;
        end
        A1 = A1 + sum_psi_tjk * spdiags(OutEff(:,t) .* beta5(:,t), 0, ns, ns) * Litj0;
        A2 = A2 + sum_psi_tjk * (spdiags(OutEff(:,t) .* beta5(:,t), 0, ns, ns) * Mgtj0 + spdiags(OutEff(:,t) .* beta4(:,t), 0, ns, ns) * Gtj0);
        A3 = A3 + sum_psi_tjk * (spdiags(OutEff(:,t) .* beta5(:,t), 0, ns, ns) * Mhtj0 + spdiags(OutEff(:,t) .* beta4(:,t), 0, ns, ns) * Htj0);
        for k = 2:Istr.idx.nc(t,j)+1
          ii  = (1:ns)';
          jj1 = (vv.i1.Psc(t,j,k):vv.iN.Psc(t,j,k))';
          jj2 = (vv.i1.Psd(t,j,k):vv.iN.Psd(t,j,k))';
          Gtjk = sparse(ii, jj1, -Istr.Delta_T  *  InEff(:,t), ns, nvars);
          Htjk = sparse(ii, jj2, -Istr.Delta_T ./ OutEff(:,t), ns, nvars);
          A6 = A6 + Istr.CostWeights(k,j,t) * spdiags(OutEff(:,t) .* beta3(:,t), 0, ns, ns) * Gtjk;
          A7 = A7 + Istr.CostWeights(k,j,t) * spdiags(OutEff(:,t) .* beta3(:,t), 0, ns, ns) * Htjk;
        end
      end
    end
    Cfstor = -baseMVA * ...
       (Istr.Storage.TerminalStoragePrice'      * (A2 + A3) + ...
        Istr.Storage.TerminalChargingPrice0'    * A4 + ...
        Istr.Storage.TerminalDischargingPrice0' * A5 + ...
        Istr.Storage.TerminalChargingPriceK'    * A6 + ...
        Istr.Storage.TerminalDischargingPriceK' * A7);
    if Istr.Storage.ForceCyclicStorage
      % If the horizon model for the storage is cyclic and therefore s0 is a
      % variable, then that initial storage must come at a cost,
      % (InitialStorageCost) otherwise the optimizer will force the gratis
      % s0 up just to have (possibly) more storage left at the end.
      Cfstor(vv.i1.S0:vv.iN.S0) = ...
        Cfstor(vv.i1.S0:vv.iN.S0) + ...
            baseMVA * Istr.Storage.InitialStorageCost';
      % and the term in the final expected storage related to s0 is also
      % not constant, so must be included in the objective function
      Cfstor(vv.i1.S0:vv.iN.S0) = ...
          Cfstor(vv.i1.S0:vv.iN.S0) - ...
          baseMVA * Istr.Storage.TerminalStoragePrice' * A1;
    end
    cp = struct('Cw', Cfstor');
    om = add_costs(om, 'fstor', cp);

    % The following is a hack to make the storage state bounds tight;
    % assign them a very small cost
    om = add_costs(om, 'SpSmFudge', {nt});
    cp = struct('Cw', 1e-2 * [-ones(ns,1); ones(ns,1)]);
    for t = 1:nt
      vs = struct('name', {'Sm', 'Sp'}, 'idx', {{t}, {t}});
      om = add_costs(om, 'SpSmFudge', {t}, cp, vs);
    end
  else
    Cfstor = sparse(1, nvars);
  end

  % Plug into struct
  if verbose
    fprintf('- Assembling full set of costs.\n');
  end
  om = build_cost_params(om, 'force');
  cp = get_cost_params(om);
  Istr.QP.Cfstor = Cfstor;
  Istr.QP.H1 = cp.N' * cp.H * cp.N;
  Istr.QP.C1 = cp.N' * cp.Cw;
  Istr.QP.c1 = c1;
end     % if mpopt.most.build_model

% With all pieces of the cost in place, can proceed to build the total
% cost now.
Istr.QP.H = Istr.QP.H1;
Istr.QP.C = Istr.QP.C1;
Istr.QP.c = Istr.QP.c1;
if isfield(Istr, 'CoordCost') && ...
        (~isempty(Istr.CoordCost.Cuser) || ~isempty(Istr.CoordCost.Huser))
  if verbose
    fprintf('- Adding coordination cost to standard cost.\n');
  end
  nvuser = length(Istr.CoordCost.Cuser);
  nvars = Istr.idx.nvars;
  Istr.QP.H = Istr.QP.H + ...
            [ Istr.CoordCost.Huser       sparse(nvuser,nvars-nvuser) ;
            sparse(nvars-nvuser,nvuser)  sparse(nvars-nvuser,nvars-nvuser) ];
  Istr.QP.C(1:nvuser) = Istr.QP.C(1:nvuser) +  Istr.CoordCost.Cuser(:);
  Istr.QP.c = Istr.QP.c + Istr.CoordCost.cuser;
  
%   cp = struct('Cw', Istr.CoordCost.Cuser(:), ...
%         'H', [ Istr.CoordCost.Huser     sparse(nvuser,nvars-nvuser) ;
%             sparse(nvars-nvuser,nvuser) sparse(nvars-nvuser,nvars-nvuser) ]);
%   om = add_costs(om, 'CoordCost', cp);
%   om = build_cost_params(om, 'force');
end

Istr.om = om;
[vv, ll] = get_idx(om);
if verbose
  fprintf('- Assembling full set of constraints.\n');
end
[Istr.QP.A, Istr.QP.l, Istr.QP.u] = linear_constraints(om);
if verbose
  fprintf('- Assembling full set of variable bounds.\n');
end
[Istr.QP.x0, Istr.QP.xmin, Istr.QP.xmax, Istr.QP.vtype] = getv(om);

% cp = get_cost_params(om);
% Istr.QP.H = cp.N' * cp.H * cp.N;
% Istr.QP.C = cp.N' * cp.Cw;
% oldidx(Istr, Istr);
% Istrtmp = oldidx(Istr);
% oldidx(Istr, Istrtmp);

tmptime(2,:) = clock;

% TEST SECTION
% create generation excess variables to impose a cost on excess generation
% nexcess = Istr.idx.nf_total;
% kk = nvars + 1;
% for t = 1:nt
%   for j = 1:Istr.idx.nj(t)
%     for k = 1:Istr.idx.nc(t,j)+1
%       Istr.idx.exbase(t,j,k) = kk;
%       Istr.idx.exend(t,j,k) = kk;
%       kk = kk + 1;
%     end
%   end
% end
% nvars = kk - 1;
% Aex = sparse(0, nvars);
% lex = [];
% uex = [];
% kk = size(Istr.QP.A, 1) + 1;
% for t = 1:nt
%   for j = 1:Istr.idx.nj(t)
%     for k = 1:Istr.idx.nc(t,j)+1
%       Aex = [ Aex;
%               sparse( ones(ng+1,1), [Istr.idx.pbas(t,j,k):Istr.idx.pend(t,j,k) Istr.idx.exbase(t,j,k)]', [ones(ng,1); -1], 1, nvars) ];
%       lex = [ lex; 1.01*sum(Istr.flow(t,j,k).mpc.bus(:,PD))/baseMVA ];
%       uex = [ uex; Inf ];
%     end
%   end
% end
% Istr.QP.A = [ Istr.QP.A sparse(size(Istr.QP.A,1), nvars-Istr.idx.nvars);
%               Aex];
% Istr.QP.l = [ Istr.QP.l; lex];
% Istr.QP.u = [ Istr.QP.u; uex];
% Istr.QP.xmin = [ Istr.QP.xmin; zeros(nexcess,1) ];
% Istr.QP.xmax = [ Istr.QP.xmax; ones(nexcess,1) ];
% Istr.QP.H = [ Istr.QP.H sparse(Istr.idx.nvars, nvars-Istr.idx.nvars);
%               sparse(nvars-Istr.idx.nvars, nvars) ];
% Istr.QP.C = [ Istr.QP.C;  1e5*ones(nexcess,1)];


% Call solver!
Ostr = Istr;
if mpopt.most.solve_model
  %% check consistency of model options (in case Istr was built in previous call)
  if Istr.DCMODEL ~= mo.DCMODEL
    error('MD.DCMODEL inconsistent with MPOPT.most.dc_model');
  end
  if Istr.SecurityConstrained ~= mo.SecurityConstrained
    error('MD.SecurityConstrained inconsistent with MPOPT.most.security_constraints (and possible presence of MD.cont(t,j).contab)');
  end
  if Istr.QCoordination ~= mo.QCoordination
    error('MD.QCoordination inconsistent with MPOPT.most.q_coordination');
  end
  if Istr.Storage.ForceCyclicStorage ~= mo.ForceCyclicStorage
    error('MD.Storage.ForceCyclicStorage inconsistent with MPOPT.most.storage.cyclic');
  end
  if Istr.Storage.ForceExpectedTerminalStorage ~= mo.ForceExpectedTerminalStorage
    error('MD.Storage.ForceExpectedTerminalStorage inconsistent with MPOPT.most.storage.terminal_target (and possible presence of MD.Storage.ExpectedTerminalStorageAim|Min|Max)');
  end
  if Istr.UC.run ~= UC
    error('MD.UC.run inconsistent with MPOPT.most.uc.run (and possible presence of MD.UC.CommitKey)');
  end
  %% set options
  if any(any(Istr.QP.H))
    model = 'QP';
  else
    model = 'LP';
  end
  if UC
    model = ['MI' model];
  end
  Ostr.QP.opt = mpopt2qpopt(mpopt, model, 'most');
  if verbose
    fprintf('- Calling %s solver.\n\n', model);
    fprintf('============================================================================\n\n');
  end
  if UC
    [Ostr.QP.x, Ostr.QP.f, Ostr.QP.exitflag, Ostr.QP.output, ...
            Ostr.QP.lambda ] = miqps_matpower( Istr.QP.H, Istr.QP.C, ...
                Istr.QP.A, Istr.QP.l, Istr.QP.u, Istr.QP.xmin, Istr.QP.xmax, ...
                [], Istr.QP.vtype, Ostr.QP.opt);
  else
    [Ostr.QP.x, Ostr.QP.f, Ostr.QP.exitflag, Ostr.QP.output, ...
            Ostr.QP.lambda ] = qps_matpower( Istr.QP.H, Istr.QP.C, ...
                Istr.QP.A, Istr.QP.l, Istr.QP.u, Istr.QP.xmin, Istr.QP.xmax, ...
                [], Ostr.QP.opt);
  end
  if Ostr.QP.exitflag > 0
    if verbose
      fprintf('\n============================================================================\n');
      fprintf('- MOST: %s solved successfully.\n', model);
    end
  else
    fprintf('\n============================================================================\n');
    fprintf('- MOST: %s solver ''%s'' failed with exit flag = %d\n', model, Ostr.QP.opt.alg, Ostr.QP.exitflag);
    fprintf('  You can query the workspace to debug.\n')
    fprintf('  When finished, type the word "return" to continue.\n\n');
    keyboard;
  end
  % Unpack results
  if verbose
    fprintf('- Post-processing results.\n');
  end
  for t = 1:nt
    if UC
      Ostr.UC.CommitSched(:, t) = Ostr.QP.x(vv.i1.u(t):vv.iN.u(t));
    end
    for j = 1:Istr.idx.nj(t)
      for k = 1:Istr.idx.nc(t,j)+1
        mpc = Ostr.flow(t,j,k).mpc;     %% pull mpc from output struct
        % Some initialization of data
        if Ostr.DCMODEL
          mpc.bus(:, VM) = 1;
        end
        % Injections and shadow prices
        mpc.gen(:, PG) = baseMVA * Ostr.QP.x(vv.i1.Pg(t,j,k):vv.iN.Pg(t,j,k));
        %% need to update Qg for loads consistent w/constant power factor
        Pmin = mpc.gen(:, PMIN);
        Qmin = mpc.gen(:, QMIN);
        Qmax = mpc.gen(:, QMAX);
        ivl = find( isload(mpc.gen) & (Qmin ~= 0 | Qmax ~= 0) );
        Qlim = (Qmin(ivl) == 0) .* Qmax(ivl) + (Qmax(ivl) == 0) .* Qmin(ivl);
        mpc.gen(ivl, QG) = mpc.gen(ivl, PG) .* Qlim ./ Pmin(ivl);
        if Ostr.DCMODEL
          %% bus angles
          mpc.bus(:, VA) = (180/pi) * Ostr.QP.x(vv.i1.Va(t,j,k):vv.iN.Va(t,j,k));
          
          %% nodal prices
          price = (Ostr.QP.lambda.mu_u(ll.i1.Pmis(t,j,k):ll.iN.Pmis(t,j,k))-Ostr.QP.lambda.mu_l(ll.i1.Pmis(t,j,k):ll.iN.Pmis(t,j,k))) / baseMVA;
          mpc.bus(:, LAM_P) = price;
          
          %% line flows and line limit shadow prices
          mpc.branch(:, PF) = 0;
          mpc.branch(:, QF) = 0;
          mpc.branch(:, PT) = 0;
          mpc.branch(:, QT) = 0;
          mpc.branch(:, MU_SF) = 0;
          mpc.branch(:, MU_ST) = 0;
          ion = find(mpc.branch(:, BR_STATUS));
          %ioff = find(~mpc.branch(:, BR_STATUS));
          rows = ll.i1.Pf(t,j,k):ll.iN.Pf(t,j,k);
          cols = vv.i1.Va(t,j,k):vv.iN.Va(t,j,k);
          lf = baseMVA * (Ostr.QP.A(rows,cols) * Ostr.QP.x(cols) + Ostr.flow(t,j,k).PLsh);
          mpc.branch(ion, PF) = lf;
          mpc.branch(ion, PT) = -lf;
          mpc.branch(ion, MU_SF) = Ostr.QP.lambda.mu_u(rows) / baseMVA;
          mpc.branch(ion, MU_ST) = Ostr.QP.lambda.mu_l(rows) / baseMVA;
        else
          %% system price
          price = (Ostr.QP.lambda.mu_l(ll.i1.Pmis(t,j,k):ll.iN.Pmis(t,j,k))-Ostr.QP.lambda.mu_u(ll.i1.Pmis(t,j,k):ll.iN.Pmis(t,j,k))) / baseMVA;
          mpc.bus(:, LAM_P) = price;
        end
        if UC
          % igenon does not contain gens ousted because of a contingency or
          % a forced-off UC.CommitKey
          igenon = find(mpc.gen(:, GEN_STATUS));
          u = Ostr.QP.x(vv.i1.u(t):vv.iN.u(t));
          mpc.gen(igenon, GEN_STATUS) = u(igenon);
          gs = mpc.gen(igenon, GEN_STATUS) > 0; % gen status
          mpc.gen(:, MU_PMAX) = 0;
          mpc.gen(:, MU_PMIN) = 0;
          mpc.gen(igenon, MU_PMAX) = gs .* ...
                  Ostr.QP.lambda.mu_u(ll.i1.uPmax(t,j,k):ll.iN.uPmax(t,j,k)) / baseMVA;
          mpc.gen(igenon, MU_PMIN) = gs .* ...
                  Ostr.QP.lambda.mu_u(ll.i1.uPmin(t,j,k):ll.iN.uPmin(t,j,k)) / baseMVA;
          if Ostr.QCoordination
            mpc.gen(:, MU_QMAX) = 0;
            mpc.gen(:, MU_QMIN) = 0;
            mpc.gen(igenon, MU_QMAX) = gs .* ...
                    Ostr.QP.lambda.mu_u(ll.i1.uQmax(t,j,k):ll.iN.uQmax(t,j,k)) / baseMVA;
            mpc.gen(igenon, MU_QMIN) = gs .* ...
                    Ostr.QP.lambda.mu_u(ll.i1.uQmin(t,j,k):ll.iN.uQmin(t,j,k)) / baseMVA;
          end
        else
          gs = mpc.gen(:, GEN_STATUS) > 0;      % gen status
          mpc.gen(:, MU_PMAX) = gs .* ...
                  Ostr.QP.lambda.upper(vv.i1.Pg(t,j,k):vv.iN.Pg(t,j,k)) / baseMVA;
          mpc.gen(:, MU_PMIN) = gs .* ...
                  Ostr.QP.lambda.lower(vv.i1.Pg(t,j,k):vv.iN.Pg(t,j,k)) / baseMVA;
          if Ostr.QCoordination
            mpc.gen(:, MU_QMAX) = gs .* ...
                    Ostr.QP.lambda.upper(vv.i1.Qg(t,j,k):vv.iN.Qg(t,j,k)) / baseMVA;
            mpc.gen(:, MU_QMIN) = gs .* ...
                    Ostr.QP.lambda.lower(vv.i1.Qg(t,j,k):vv.iN.Qg(t,j,k)) / baseMVA;
          end
        end
        Ostr.flow(t,j,k).mpc = mpc;     %% stash modified mpc in output struct
        if Istr.IncludeFixedReserves
          z = zeros(ng, 1);
          Ostr.FixedReserves(t,j,k).R   = z;
          Ostr.FixedReserves(t,j,k).prc = z;
          Ostr.FixedReserves(t,j,k).mu = struct('l', z, 'u', z, 'Pmax', z);
          Ostr.FixedReserves(t,j,k).totalcost = compute_cost(om, Ostr.QP.x, 'Rcost', {t,j,k});
          r = Ostr.FixedReserves(t,j,k);
          r.R(r.igr) = Ostr.QP.x(vv.i1.R(t,j,k):vv.iN.R(t,j,k)) * baseMVA;
          for gg = r.igr
            iz = find(r.zones(:, gg));
            kk = ll.i1.Rreq(t,j,k):ll.iN.Rreq(t,j,k);
            r.prc(gg) = sum(Ostr.QP.lambda.mu_l(kk(iz))) / baseMVA;
          end
          r.mu.l(r.igr)    = Ostr.QP.lambda.lower(vv.i1.R(t,j,k):vv.iN.R(t,j,k)) / baseMVA;
          r.mu.u(r.igr)    = Ostr.QP.lambda.upper(vv.i1.R(t,j,k):vv.iN.R(t,j,k)) / baseMVA;
          r.mu.Pmax(r.igr) = Ostr.QP.lambda.mu_u(ll.i1.Pg_plus_R(t,j,k):ll.iN.Pg_plus_R(t,j,k)) / baseMVA;
          Ostr.FixedReserves(t,j,k) = r;
        end
      end
    end
    % Contract, contingency reserves, energy limits
    Ostr.results.Pc(:,t)  = baseMVA * Ostr.QP.x(vv.i1.Pc(t):vv.iN.Pc(t));
    Ostr.results.Rpp(:,t) = baseMVA * Ostr.QP.x(vv.i1.Rpp(t):vv.iN.Rpp(t));
    Ostr.results.Rpm(:,t) = baseMVA * Ostr.QP.x(vv.i1.Rpm(t):vv.iN.Rpm(t));
    if ns
      Ostr.results.Sm(:,t)  = baseMVA * Ostr.QP.x(vv.i1.Sm(t):vv.iN.Sm(t));
      Ostr.results.Sp(:,t)  = baseMVA * Ostr.QP.x(vv.i1.Sp(t):vv.iN.Sp(t));
    end
  end
  % Ramping reserves
  for t = 1:Ostr.idx.ntramp
    Ostr.results.Rrp(:,t) = baseMVA * Ostr.QP.x(vv.i1.Rrp(t):vv.iN.Rrp(t));
    Ostr.results.Rrm(:,t) = baseMVA * Ostr.QP.x(vv.i1.Rrm(t):vv.iN.Rrm(t));
  end
  % Expected energy prices for generators, per generator and per period,
  % both absolute and conditional on making it to that period
  Ostr.results.GenPrices = zeros(ng, nt);
  Ostr.results.CondGenPrices = zeros(ng, nt);
  for t = 1:nt
    pp = zeros(ng,1);
    for j = 1:Ostr.idx.nj(t)
      for k = 1:Ostr.idx.nc(t,j)+1
        pp = pp + Ostr.flow(t,j,k).mpc.bus(Ostr.flow(t,j,k).mpc.gen(:,GEN_BUS), LAM_P);
      end
    end
    Ostr.results.GenPrices(:,t) = pp;
    Ostr.results.CondGenPrices(:, t) = pp / Ostr.StepProb(t);
  end
  % Obtain contingency reserve prices, per generator and period
  Ostr.results.RppPrices = zeros(ng, nt);
  Ostr.results.RpmPrices = zeros(ng, nt);
  for t = 1:nt
    Ostr.results.RppPrices(:, t) = Ostr.QP.lambda.lower(vv.i1.Rpp(t):vv.iN.Rpp(t)) / baseMVA;
    Ostr.results.RpmPrices(:, t) = Ostr.QP.lambda.lower(vv.i1.Rpm(t):vv.iN.Rpm(t)) / baseMVA;
    for j = 1:Istr.idx.nj(t);
      for k = 1:Istr.idx.nc(t,j)+1
        ii = find(Istr.flow(t,j,k).mpc.gen(:, GEN_STATUS) > 0);
        Ostr.results.RppPrices(ii, t) = Ostr.results.RppPrices(ii, t) + Ostr.QP.lambda.mu_l(ll.i1.dPpRp(t,j,k):ll.iN.dPpRp(t,j,k)) / baseMVA;
        Ostr.results.RpmPrices(ii, t) = Ostr.results.RpmPrices(ii, t) + Ostr.QP.lambda.mu_l(ll.i1.dPmRm(t,j,k):ll.iN.dPmRm(t,j,k)) / baseMVA;
      end
    end
  end
  % Obtain ramping reserve prices, per generator and period
  Ostr.results.RrpPrices = zeros(ng, Ostr.idx.ntramp);
  Ostr.results.RrmPrices = zeros(ng, Ostr.idx.ntramp);
  % First, 1:nt-1
  for t = 1:nt-1
    for j1 = 1:Ostr.idx.nj(t)
      for j2 = 1:Ostr.idx.nj(t+1)
        if Istr.tstep(t+1).TransMask(j2,j1)
          Ostr.results.RrpPrices(:, t) = Ostr.results.RrpPrices(:, t) + Ostr.QP.lambda.mu_l(ll.i1.Rrp(t,j1,j2):ll.iN.Rrp(t,j1,j2)) / baseMVA;
          Ostr.results.RrmPrices(:, t) = Ostr.results.RrmPrices(:, t) + Ostr.QP.lambda.mu_l(ll.i1.Rrm(t,j1,j2):ll.iN.Rrm(t,j1,j2)) / baseMVA;
        end
      end
    end
  end
  % then last period only if specified for with terminal state
  if ~Ostr.OpenEnded
    for j1 = 1:Ostr.idx.nj(nt)
      Ostr.results.RrpPrices(:, nt) = Ostr.results.RrpPrices(:, nt) + Ostr.QP.lambda.mu_l(ll.i1.Rrp(nt,j1,1):ll.iN.Rrp(nt,j1,1)) / baseMVA;
      Ostr.results.RrmPrices(:, nt) = Ostr.results.RrmPrices(:, nt) + Ostr.QP.lambda.mu_l(ll.i1.Rrm(nt,j1,1):ll.iN.Rrm(nt,j1,1)) / baseMVA;
    end
  end
  % Expected wear and tear costs per gen and period
  Ostr.results.ExpectedRampCost = zeros(ng, Ostr.idx.ntramp+1);
  % First do first period wrt to InitialPg.
  for j = 1:Istr.idx.nj(1)
    w = Ostr.tstep(1).TransMat(j,1);  % the probability of going from initial state to jth
    Ostr.results.ExpectedRampCost(:, 1) = Ostr.results.ExpectedRampCost(:, 1) ...
        + 0.5 * w * Ostr.RampWearCostCoeff(:,1) .* (Ostr.flow(1,j,1).mpc.gen(:,PG) - Ostr.InitialPg).^2;
  end
  % Then the remaining periods
  for t = 2:nt
    for j2 = 1:Ostr.idx.nj(t)
      for j1 = 1:Ostr.idx.nj(t-1)
        w = Ostr.tstep(t).TransMat(j2,j1) * Ostr.CostWeights(1, j1, t-1);
        Ostr.results.ExpectedRampCost(:, t) = Ostr.results.ExpectedRampCost(:, t) ...
            + 0.5 * w * Ostr.RampWearCostCoeff(:,t) .* (Ostr.flow(t,j2,1).mpc.gen(:,PG) - Ostr.flow(t-1,j1,1).mpc.gen(:,PG)) .^2;
      end
    end
  end
  % Finally, if there is a terminal state problem, apply cost to
  if ~Ostr.OpenEnded
    for j = 1:Istr.idx.nj(nt)
      w = Istr.tstep(t+1).TransMat(1, j) * Istr.CostWeights(1, j, nt);
      Ostr.results.ExpectedRampCost(:, nt+1) = 0.5 * w * Ostr.RampWearCostCoeff(:,nt+1) .* (Ostr.TerminalPg - Ostr.flow(nt,j,1).mpc.gen(:,PG)) .^2;
    end
  end
  % Compute expected dispatch, conditional on making it to the
  % corresponding period
  Ostr.results.ExpectedDispatch = zeros(ng, nt);
  for t = 1:nt
    pp = sum(Ostr.CostWeights(1,1:Ostr.idx.nj(t),t)');  % gamma(t+1)
    for j = 1:Ostr.idx.nj(t)
      Ostr.results.ExpectedDispatch(:,t) = Ostr.results.ExpectedDispatch(:,t) + ...
            Ostr.CostWeights(1,j,t)/pp * Ostr.flow(t,j,1).mpc.gen(:,PG);
    end
  end
  % If Cyclic storage, pull InitialStorage value out of x
  if ns && Ostr.Storage.ForceCyclicStorage
    Ostr.Storage.InitialStorage = baseMVA * Ostr.QP.x(vv.i1.S0:vv.iN.S0);
  end
  % Compute expected storage state trajectory
  Ostr.Storage.ExpectedStorageState = zeros(ns,nt);
  if ns
    for t = 1:nt
      pp = sum(Ostr.CostWeights(1,1:Ostr.idx.nj(t),t)');    %% gamma(t+1)
      for j = 1:Ostr.idx.nj(t)
        Lfj = Ostr.tstep(t).Lf( j:Ostr.idx.nj(t):(ns-1)*Ostr.idx.nj(t)+j, :);
        Ngj = Ostr.tstep(t).Ng( j:Ostr.idx.nj(t):(ns-1)*Ostr.idx.nj(t)+j, :);
        Nhj = Ostr.tstep(t).Nh( j:Ostr.idx.nj(t):(ns-1)*Ostr.idx.nj(t)+j, :);
        Ostr.Storage.ExpectedStorageState(:,t) = ...
            Ostr.Storage.ExpectedStorageState(:,t) + ...
                baseMVA * Ostr.CostWeights(1,j,t)/pp * ...
                ( Lfj * Ostr.Storage.InitialStorage/baseMVA + ...
                  (Ngj + Nhj) * Ostr.QP.x );
      end
    end
    Ostr.Storage.ExpectedStorageDispatch = ...
        Ostr.results.ExpectedDispatch(Ostr.Storage.UnitIdx, :);
  end
  % If there is a dynamical system, extract the state vectors and outputs
  % from the solution
  if nyt
    if nzd
      Ostr.results.Z = zeros(nzd, nyt);
      for t = 1:nyt
        Ostr.results.Z(:,t) = Ostr.QP.x(vv.i1.Z(t):vv.iN.Z(t));
      end
    end
    Ostr.results.Y = zeros(nyo, nyt);
    if nyo
      for t = 1:nyt
        Ostr.results.Y(:, t) = ...
                Ostr.QP.A(ll.i1.DSy(t):ll.iN.DSy(t), :) * Ostr.QP.x;
      end
    end
  end
  Ostr.results.f = Ostr.QP.f;
end % if mpopt.most.solve_model

tmptime(3,:) = clock;

Ostr.results.SetupTime = etime(tmptime(2,:), tmptime(1,:));
Ostr.results.SolveTime = etime(tmptime(3,:), tmptime(2,:));

if verbose
  fprintf('- MOST: Done.\n\n');
end

return;

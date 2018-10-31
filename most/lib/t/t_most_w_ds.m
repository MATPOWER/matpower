function t_most_w_ds(quiet)
%T_MOST_W_DS  Test for MOST with dynamical system constraints.

%   MOST
%   Copyright (c) 2015-2016, Power Systems Engineering Research Center (PSERC)
%   by Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Nacional de Colombia
%   and Ray Zimmerman, PSERC Cornell
%
%   This file is part of MOST.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/most for more info.

if nargin < 1
    quiet = 0;
end

include_MIPS = 0;   %% set to 1, to attempt even if MIPS is the best solver
                    %% available (takes a LONG time and currently fails)
n_tests = 1;

t_begin(n_tests, quiet);

casefile = 'c118swf';
solnfile = 't_most_w_ds_z';

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

if have_fcn('cplex') || have_fcn('gurobi') || have_fcn('mosek') || ...
        have_fcn('quadprog_ls') || include_MIPS
    mdi = md_init;

    mpopt = mpoption('verbose', 0);

    %% choose solver
    if have_fcn('cplex')
        mpopt = mpoption(mpopt, 'most.solver', 'CPLEX');
    elseif have_fcn('gurobi')
        mpopt = mpoption(mpopt, 'most.solver', 'GUROBI');
    elseif have_fcn('quadprog_ls')
        mpopt = mpoption(mpopt, 'most.solver', 'OT');
    elseif have_fcn('mosek')
        mpopt = mpoption(mpopt, 'most.solver', 'MOSEK');
    else
        mpopt = mpoption(mpopt, 'most.solver', 'MIPS');
        if have_fcn('pardiso')
            mpopt = mpoption(mpopt, 'mips.linsolver', 'PARDISO');
        end
    end

    %% set options
    if have_fcn('cplex')
        mpopt = mpoption(mpopt, 'cplex.opts.threads', 2);   % set this manually here
    end
    if have_fcn('gurobi')
        mpopt = mpoption(mpopt, 'most.solver', 'GUROBI');
        mpopt = mpoption(mpopt, 'gurobi.method', 2);        %% barrier
        mpopt = mpoption(mpopt, 'gurobi.threads', 2);
        mpopt = mpoption(mpopt, 'gurobi.opts.BarConvTol', 1e-6);        %% 1e-8
        mpopt = mpoption(mpopt, 'gurobi.opts.FeasibilityTol', 1e-4);    %% 1e-6
        mpopt = mpoption(mpopt, 'gurobi.opts.OptimalityTol', 1e-5);     %% 1e-6
    end
    if have_fcn('quadprog_ls')
        mpopt = mpoption(mpopt, 'most.solver', 'OT');
        mpopt = mpoption(mpopt, 'quadprog.TolFun', 1e-13);
    end
    if have_fcn('mosek')
        mpopt = mpoption(mpopt, 'most.solver', 'MOSEK');
        mpopt = mpoption(mpopt, 'mosek.num_threads', 2);
    else
        mpopt = mpoption(mpopt, 'most.solver', 'MIPS');
        mpopt = mpoption(mpopt, 'mips.max_it', 500);
        if have_fcn('pardiso')
            mpopt = mpoption(mpopt, 'mips.linsolver', 'PARDISO');
        end
    end
    % mpopt = mpoption(mpopt, 'most.solver', 'GUROBI');
    % mpopt = mpoption(mpopt, 'most.solver', 'CLP');
    % mpopt = mpoption(mpopt, 'most.solver', 'IPOPT');
    % mpopt = mpoption(mpopt, 'most.solver', 'MIPS');

    mdi.mpc = loadcase(casefile);
    mdi.InitialPg = mdi.mpc.gen(:,PG);
    nt = 24;
    ng = size(mdi.mpc.gen, 1);
    mdi.idx.nt = nt;
    PositiveActiveReservePrice = ones(ng,1);
    PositiveActiveReserveQuantity = 0.25*mdi.mpc.gen(:,PMAX);
    NegativeActiveReservePrice = ones(ng,1);
    NegativeActiveReserveQuantity = PositiveActiveReserveQuantity;
    PositiveActiveDeltaPrice = ones(ng,1);
    NegativeActiveDeltaPrice = ones(ng,1);
    PositiveLoadFollowReservePrice = ones(ng,1);
    PositiveLoadFollowReserveQuantity = 0.5*mdi.mpc.gen(:,PMAX);
    NegativeLoadFollowReservePrice = ones(ng,1);
    NegativeLoadFollowReserveQuantity = PositiveLoadFollowReserveQuantity;
    %mdi.mpc.gen(:,RAMP_10) = 0.20 * mdi.mpc.gen(PMAX);
    %mdi.mpc.gen(:,RAMP_AGC) = 0.20 * mdi.mpc.gen(PMAX);
    %mdi.mpc.gen(:,RAMP_30) = 0.50 * mdi.mpc.gen(PMAX);
    mdi.mpc.gen(:,RAMP_10) = 1.0 * mdi.mpc.gen(:,PMAX);
    mdi.mpc.gen(:,RAMP_AGC) = 1.0 * mdi.mpc.gen(:,PMAX);
    mdi.mpc.gen(:,RAMP_30) = 1.0 * mdi.mpc.gen(:,PMAX);

    mdi.RampWearCostCoeff = 0.05 * ones(ng,1);   % (i, t) note different scheme!
    for t = 2:nt
      mdi.RampWearCostCoeff(:, t) = mdi.RampWearCostCoeff(:, 1);
    end
    mdi.Storage(1).UnitIdx = mdi.mpc.iess;
    ns = length(mdi.Storage.UnitIdx);
    Minstor = zeros(ns,1);
    Maxstor = 200 * ones(ns,1);
    %mdi.Storage.MinStorageLevel    = zeros(ns,1);
    %mdi.Storage.MaxStorageLevel    = 200 * ones(ns,1);
    mdi.Storage.InitialStorage     = 50 * ones(ns,1);
    mdi.Storage.InitialStorageLowerBound = 50*ones(ns,1);
    mdi.Storage.InitialStorageUpperBound = 50*ones(ns,1);
    mdi.Storage.OutEff             = 0.95 * ones(ns,1);
    mdi.Storage.InEff              = 0.9  * ones(ns ,1);
    mdi.Storage.InitialStorageCost         = 35 * ones(ns, 1);
    mdi.Storage.TerminalStoragePrice       = 35 * ones(ns, 1); % applied to psc_tij0, psd_tij0 (non-terminal states)
    mdi.Storage.TerminalChargingPrice0     = 35 * ones(ns, 1); % applied to psc_tijk (contingency terminal states)
    mdi.Storage.TerminalDischargingPrice0  = 35 * ones(ns, 1); % applied to psd_tijk (contingency terminal states)
    mdi.Storage.TerminalChargingPriceK     = 10 * ones(ns, 1); % applied to psc_tij0 (end-of-horizon terminal states)
    mdi.Storage.TerminalDischargingPriceK  = 40 * ones(ns, 1); % applied to psd_tij0 (end-of-horizon terminal states)
    mpopt = mpoption(mpopt, 'most.storage.terminal_target', 0);
    mdi.Storage.ExpectedTerminalStorageAim = mdi.Storage.InitialStorage;  % expected terminal storage if mpopt.most.storage.terminal_target is true
    mdi.Storage.LossFactor         = zeros(ns,1);  % fraction of storage lost in each period
    mdi.Storage.IncludeValueOfTerminalStorage = 1;
    mpopt = mpoption(mpopt, 'most.storage.cyclic', 1);

    for t = 1:nt
      mdi.offer(t).gencost = mdi.mpc.gencost;
      mdi.offer(t).PositiveActiveReservePrice = PositiveActiveReservePrice;
      mdi.offer(t).PositiveActiveReserveQuantity = PositiveActiveReserveQuantity;
      mdi.offer(t).NegativeActiveReservePrice = NegativeActiveReservePrice;
      mdi.offer(t).NegativeActiveReserveQuantity = NegativeActiveReserveQuantity;
      mdi.offer(t).PositiveActiveDeltaPrice = PositiveActiveDeltaPrice;
      mdi.offer(t).NegativeActiveDeltaPrice = NegativeActiveDeltaPrice;
      mdi.offer(t).PositiveLoadFollowReservePrice = PositiveLoadFollowReservePrice;
      mdi.offer(t).PositiveLoadFollowReserveQuantity = PositiveLoadFollowReserveQuantity;
      mdi.offer(t).NegativeLoadFollowReservePrice = NegativeLoadFollowReservePrice;
      mdi.offer(t).NegativeLoadFollowReserveQuantity = NegativeLoadFollowReserveQuantity;
      mdi.Storage.MinStorageLevel(:,t) = Minstor;
      mdi.Storage.MaxStorageLevel(:,t) = Maxstor;
    end
    mdi.Storage.MinStorageLevel(:,nt+1) = Minstor;  % Needed if mpopt.most.storage.cyclic
    mdi.Storage.MaxStorageLevel(:,nt+1) = Maxstor;

    %mdi.Storage.MinStorageLevel(:,4) = [10 ; 10];  % is this enough to create infeasibility?
    %mdi.Storage.MaxStorageLevel(:,4) = [10; 50 ];

    mdi.UC.CommitSched = ones(ng,nt);

    mdi.Delta_T = 1;
    %            0:0  1:00 2:00 3:00 4:00  5:00  6:00  7:00 8:00 9:00 10:00 11:00 12:00 13:00 14:00 15:00 16:00 17:00 18:00 19:00 20:00 21:00 22:00 23:00 
    loadprof = [ 0.6  0.6   0.6 0.7  0.75  0.8   0.9   1.1  1.2  1.3   1.4  1.4   1.3   1.4    1.4   1.4   1.4   1.3   1.1   1.1   1.0   0.9   0.8   0.7  ];


    %                label  probability   type      row     column      chg type  newvalue
    partialcontabrow =[ 1       0        CT_TBUS     0        PD         CT_REL ];
    %mdi.tstep(1).OpCondSched(1).tab= [ 
    %                   1        0        CT_TBUS     0        PD         CT_REL     0.8  ];
    %mdi.tstep(2).OpCondSched(1).tab= [ 
    %                   1        0        CT_TBUS     0        PD         CT_REL     1.0  ];
    %mdi.tstep(3).OpCondSched(1).tab= [ 
    %                   1        0        CT_TBUS     0        PD         CT_REL     1.2  ];
    %mdi.tstep(4).OpCondSched(1).tab= [
    %                   1        0        CT_TBUS     0        PD         CT_REL     0.9 ];

    for t = 1:nt
      mdi.tstep(t).OpCondSched(1).tab = [ partialcontabrow   loadprof(t) ];
    end


    for t = 1:nt
       mdi.tstep(t).OpCondSched(2).tab = mdi.tstep(t).OpCondSched(1).tab;
       mdi.tstep(t).OpCondSched(2).tab(1,7) = 1.1*mdi.tstep(t).OpCondSched(2).tab(1,7);
       mdi.tstep(t).OpCondSched(3).tab = mdi.tstep(t).OpCondSched(1).tab;
       mdi.tstep(t).OpCondSched(3).tab(1,7) = 0.9*mdi.tstep(t).OpCondSched(1).tab(1,7);
    end


    contab = [%         1       0.01      CT_TBUS     0        PD         CT_REL     1.05 ;
                       1       0.01      CT_TGEN     2     GEN_STATUS    CT_REP      0    ;
                       2       0.01      CT_TGEN     5     GEN_STATUS    CT_REP      0    ;
                       ];

    for t = 1:nt
      for j = 1:3  % mdi.idx.nj(t)
        mdi.cont(t,j).contab = contab;
      end
    end



    mdi.tstep(1).TransMat = [ 1/3;
                               1/3
                               1/3];
    for t = 2:nt
      mdi.tstep(t).TransMat = 1/3 * ones(3,3);
    end


    ntds = 24;
    mdi.idx.ntds = ntds;
    m1 = 8;
    m2 = 12;
    B = sparse(m1*m2, ng);
    ilist = [ 2 3 4 5   3 4 5 6   3 4 5 7   3 4 5 7   3 4 5 6   4 5 6 7    2 3 4 ];
    jlist = [ 2 2 2 2   3 3 3 3   5 5 5 5   6 6 6 6   7 7 7 7   8 8 8 8    11 11 11 ];
    for i = 1:length(mdi.mpc.icoal)
     B((jlist(i)-1)*m1+ilist(i), mdi.mpc.icoal(i)) = 0.1;
    end
    A = mkdif(m1, m2, 0.5, 0.97, [1.0 0]);
    C = [];
    D = [];
    zmin = zeros(m1*m2, 1);
    zmax = 100*ones(m1*m2, 1);
    ymin = 0;
    ymax = 100;
    for t = 1:ntds
     mdi.dstep(t).A = A;
     mdi.dstep(t).B = B;
     mdi.dstep(t).C = C;
     mdi.dstep(t).D = D;
     mdi.dstep(t).zmin = zmin;
     mdi.dstep(t).zmax = zmax;
     mdi.dstep(t).ymin = ymin;
     mdi.dstep(t).ymax = ymax;
    end
    mdi.z1 = zeros(m1*m2, 1);

    mdo = most(mdi, mpopt);

    s = load(solnfile);

    t = 'dynamical system state (Z)';
    t_is(mdo.results.Z, s.Z, 4, t);
else
    t_skip(1, 'requires CPLEX, Gurobi, MOSEK or quadprog');
end

% YorN = input('Play movie? (y/n) : ', 's');
% if strcmp(upper(YorN(1)), 'Y')
%     domovie;
% end

t_end

function A = mkdif(m1, m2, alpha, r, w)
% A = mkdif(m1, m2, r, w)
%
% computes an A matrix for the difussion equations in an m1 x m2 grid,
% using a difussion speed factor alpha <= 1, a dissipation factor r < 1 and 
% a "wind" vector w = [wx wy] that roughly tells where the pollutants go.

% Carlos Murillo.  As naive as can be.

if norm(w) > 1
  w = w / norm(w);
end

n = m1 * m2;
A = sparse(n, n);
north = alpha / 4;
east = north;
south = north;
west = north;
self = (1-alpha);
if w(1) > 0
  west = west + w(1)*self;
  self = (1-w(1)) * self;
elseif w(2) < 0
  east = east + (-w(1))*self;
  self = (1+w(1)) * self;
end
if w(2) > 0
  south = south + w(2)*self;
  self = (1-w(2)) * self;
elseif w(2) < 0
  north = north + (-w(2))*self;
  self = (1+w(2)) * self;
end
tot = self + north + south + east + west;
self = (r/tot) * self;
north = (r/tot) * north;
south = (r/tot) * south;
east = (r/tot) * east;
west = (r/tot) * west;

for i = 1:m1
  for j = 1:m2
    A((j-1)*m1+i, (j-1)*m1+i) = self;
    % North
    if i > 1
      A((j-1)*m1+i, (j-1)*m1+i-1) = north;
    end
    % South
    if i < m1
      A((j-1)*m1+i, (j-1)*m1+i+1) = south;
    end
    % West
    if j > 1
      A((j-1)*m1+i, (j-2)*m1+i) = west;
    end
    % East
    if j < m2
      A((j-1)*m1+i, j*m1+i) = east;
    end
  end
end

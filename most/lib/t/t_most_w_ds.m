function t_most_w_ds(quiet)
%T_MOST_W_DS  Test for MOST with dynamical system constraints.

if nargin < 1
    quiet = 0;
end

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
        have_fcn('quadprog')
    Istr = md_init;

    mpopt = mpoption('verbose', 0);

    if have_fcn('cplex')
        mpopt = mpoption(mpopt, 'most.solver', 'CPLEX');
        mpopt = mpoption(mpopt, 'cplex.opts.threads', 2);   % set this manually here
    elseif have_fcn('gurobi')
        mpopt = mpoption(mpopt, 'most.solver', 'GUROBI');
        mpopt = mpoption(mpopt, 'gurobi.method', 2);        %% barrier
        mpopt = mpoption(mpopt, 'gurobi.threads', 2);
        mpopt = mpoption(mpopt, 'gurobi.opts.BarConvTol', 1e-7);        %% 1e-8
        mpopt = mpoption(mpopt, 'gurobi.opts.FeasibilityTol', 1e-5);    %% 1e-6
        mpopt = mpoption(mpopt, 'gurobi.opts.OptimalityTol', 1e-5);     %% 1e-6
    elseif have_fcn('quadprog')
        mpopt = mpoption(mpopt, 'most.solver', 'OT');
        mpopt = mpoption(mpopt, 'quadprog.TolFun', 1e-13);
    elseif have_fcn('mosek')
        mpopt = mpoption(mpopt, 'most.solver', 'MOSEK');
        mpopt = mpoption(mpopt, 'mosek.num_threads', 2);
    end
    % mpopt = mpoption(mpopt, 'most.solver', 'CLP');
    % mpopt = mpoption(mpopt, 'most.solver', 'IPOPT');
    % mpopt = mpoption(mpopt, 'most.solver', 'MIPS');
    % mpopt = mpoption(mpopt, 'mips.linsolver', 'PARDISO');

    Istr.mpc = loadcase(casefile);
    Istr.InitialPg = Istr.mpc.gen(:,PG);
    nt = 24;
    ng = size(Istr.mpc.gen, 1);
    Istr.idx.nt = nt;
    PositiveActiveReservePrice = ones(ng,1);
    PositiveActiveReserveQuantity = 0.25*Istr.mpc.gen(:,PMAX);
    NegativeActiveReservePrice = ones(ng,1);
    NegativeActiveReserveQuantity = PositiveActiveReserveQuantity;
    PositiveActiveDeltaPrice = ones(ng,1);
    NegativeActiveDeltaPrice = ones(ng,1);
    PositiveLoadFollowReservePrice = ones(ng,1);
    PositiveLoadFollowReserveQuantity = 0.5*Istr.mpc.gen(:,PMAX);
    NegativeLoadFollowReservePrice = ones(ng,1);
    NegativeLoadFollowReserveQuantity = PositiveLoadFollowReserveQuantity;
    %Istr.mpc.gen(:,RAMP_10) = 0.20 * Istr.mpc.gen(PMAX);
    %Istr.mpc.gen(:,RAMP_AGC) = 0.20 * Istr.mpc.gen(PMAX);
    %Istr.mpc.gen(:,RAMP_30) = 0.50 * Istr.mpc.gen(PMAX);
    Istr.mpc.gen(:,RAMP_10) = 1.0 * Istr.mpc.gen(:,PMAX);
    Istr.mpc.gen(:,RAMP_AGC) = 1.0 * Istr.mpc.gen(:,PMAX);
    Istr.mpc.gen(:,RAMP_30) = 1.0 * Istr.mpc.gen(:,PMAX);

    Istr.RampWearCostCoeff = 0.05 * ones(ng,1);   % (i, t) note different scheme!
    for t = 2:nt
      Istr.RampWearCostCoeff(:, t) = Istr.RampWearCostCoeff(:, 1);
    end
    Istr.Storage(1).UnitIdx = Istr.mpc.iess;
    ns = length(Istr.Storage.UnitIdx);
    Minstor = zeros(ns,1);
    Maxstor = 200 * ones(ns,1);
    %Istr.Storage.MinStorageLevel    = zeros(ns,1);
    %Istr.Storage.MaxStorageLevel    = 200 * ones(ns,1);
    Istr.Storage.InitialStorage     = 50 * ones(ns,1);
    Istr.Storage.InitialStorageLowerBound = 50*ones(ns,1);
    Istr.Storage.InitialStorageUpperBound = 50*ones(ns,1);
    Istr.Storage.OutEff             = 0.95 * ones(ns,1);
    Istr.Storage.InEff              = 0.9  * ones(ns ,1);
    Istr.Storage.InitialStorageCost         = 35 * ones(ns, 1);
    Istr.Storage.TerminalStoragePrice       = 35 * ones(ns, 1); % applied to psc_tij0, psd_tij0 (non-terminal states)
    Istr.Storage.TerminalChargingPrice0     = 35 * ones(ns, 1); % applied to psc_tijk (contingency terminal states)
    Istr.Storage.TerminalDischargingPrice0  = 35 * ones(ns, 1); % applied to psd_tijk (contingency terminal states)
    Istr.Storage.TerminalChargingPriceK     = 10 * ones(ns, 1); % applied to psc_tij0 (end-of-horizon terminal states)
    Istr.Storage.TerminalDischargingPriceK  = 40 * ones(ns, 1); % applied to psd_tij0 (end-of-horizon terminal states)
    mpopt = mpoption(mpopt, 'most.storage.terminal_target', 0);
    Istr.Storage.ExpectedTerminalStorageAim = Istr.Storage.InitialStorage;  % expected terminal storage if mpopt.most.storage.terminal_target is true
    Istr.Storage.LossFactor         = zeros(ns,1);  % fraction of storage lost in each period
    Istr.Storage.IncludeValueOfTerminalStorage = 1;
    mpopt = mpoption(mpopt, 'most.storage.cyclic', 1);

    for t = 1:nt
      Istr.offer(t).gencost = Istr.mpc.gencost;
      Istr.offer(t).PositiveActiveReservePrice = PositiveActiveReservePrice;
      Istr.offer(t).PositiveActiveReserveQuantity = PositiveActiveReserveQuantity;
      Istr.offer(t).NegativeActiveReservePrice = NegativeActiveReservePrice;
      Istr.offer(t).NegativeActiveReserveQuantity = NegativeActiveReserveQuantity;
      Istr.offer(t).PositiveActiveDeltaPrice = PositiveActiveDeltaPrice;
      Istr.offer(t).NegativeActiveDeltaPrice = NegativeActiveDeltaPrice;
      Istr.offer(t).PositiveLoadFollowReservePrice = PositiveLoadFollowReservePrice;
      Istr.offer(t).PositiveLoadFollowReserveQuantity = PositiveLoadFollowReserveQuantity;
      Istr.offer(t).NegativeLoadFollowReservePrice = NegativeLoadFollowReservePrice;
      Istr.offer(t).NegativeLoadFollowReserveQuantity = NegativeLoadFollowReserveQuantity;
      Istr.Storage.MinStorageLevel(:,t) = Minstor;
      Istr.Storage.MaxStorageLevel(:,t) = Maxstor;
    end
    Istr.Storage.MinStorageLevel(:,nt+1) = Minstor;  % Needed if mpopt.most.storage.cyclic
    Istr.Storage.MaxStorageLevel(:,nt+1) = Maxstor;

    %Istr.Storage.MinStorageLevel(:,4) = [10 ; 10];  % is this enough to create infeasibility?
    %Istr.Storage.MaxStorageLevel(:,4) = [10; 50 ];

    Istr.UC.CommitSched = ones(ng,nt);

    Istr.Delta_T = 1;
    %            0:0  1:00 2:00 3:00 4:00  5:00  6:00  7:00 8:00 9:00 10:00 11:00 12:00 13:00 14:00 15:00 16:00 17:00 18:00 19:00 20:00 21:00 22:00 23:00 
    loadprof = [ 0.6  0.6   0.6 0.7  0.75  0.8   0.9   1.1  1.2  1.3   1.4  1.4   1.3   1.4    1.4   1.4   1.4   1.3   1.1   1.1   1.0   0.9   0.8   0.7  ];


    %                label  probability   type      row     column      chg type  newvalue
    partialcontabrow =[ 1       0        CT_TBUS     0        PD         CT_REL ];
    %Istr.tstep(1).OpCondSched(1).tab= [ 
    %                   1        0        CT_TBUS     0        PD         CT_REL     0.8  ];
    %Istr.tstep(2).OpCondSched(1).tab= [ 
    %                   1        0        CT_TBUS     0        PD         CT_REL     1.0  ];
    %Istr.tstep(3).OpCondSched(1).tab= [ 
    %                   1        0        CT_TBUS     0        PD         CT_REL     1.2  ];
    %Istr.tstep(4).OpCondSched(1).tab= [
    %                   1        0        CT_TBUS     0        PD         CT_REL     0.9 ];

    for t = 1:nt
      Istr.tstep(t).OpCondSched(1).tab = [ partialcontabrow   loadprof(t) ];
    end


    for t = 1:nt
       Istr.tstep(t).OpCondSched(2).tab = Istr.tstep(t).OpCondSched(1).tab;
       Istr.tstep(t).OpCondSched(2).tab(1,7) = 1.1*Istr.tstep(t).OpCondSched(2).tab(1,7);
       Istr.tstep(t).OpCondSched(3).tab = Istr.tstep(t).OpCondSched(1).tab;
       Istr.tstep(t).OpCondSched(3).tab(1,7) = 0.9*Istr.tstep(t).OpCondSched(1).tab(1,7);
    end


    contab = [%         1       0.01      CT_TBUS     0        PD         CT_REL     1.05 ;
                       1       0.01      CT_TGEN     2     GEN_STATUS    CT_REP      0    ;
                       2       0.01      CT_TGEN     5     GEN_STATUS    CT_REP      0    ;
                       ];

    for t = 1:nt
      for j = 1:3  % Istr.idx.nj(t)
        Istr.cont(t,j).contab = contab;
      end
    end



    Istr.tstep(1).TransMat = [ 1/3;
                               1/3
                               1/3];
    for t = 2:nt
      Istr.tstep(t).TransMat = 1/3 * ones(3,3);
    end


    nyt = 24;
    Istr.idx.nyt = nyt;
    m1 = 8;
    m2 = 12;
    B = sparse(m1*m2, ng);
    ilist = [ 2 3 4 5   3 4 5 6   3 4 5 7   3 4 5 7   3 4 5 6   4 5 6 7    2 3 4 ];
    jlist = [ 2 2 2 2   3 3 3 3   5 5 5 5   6 6 6 6   7 7 7 7   8 8 8 8    11 11 11 ];
    for i = 1:length(Istr.mpc.icoal)
     B((jlist(i)-1)*m1+ilist(i), Istr.mpc.icoal(i)) = 0.1;
    end
    A = mkdif(m1, m2, 0.5, 0.97, [1.0 0]);
    C = [];
    D = [];
    zmin = zeros(m1*m2, 1);
    zmax = 100*ones(m1*m2, 1);
    ymin = 0;
    ymax = 100;
    for t = 1:nyt
     Istr.dstep(t).A = A;
     Istr.dstep(t).B = B;
     Istr.dstep(t).C = C;
     Istr.dstep(t).D = D;
     Istr.dstep(t).zmin = zmin;
     Istr.dstep(t).zmax = zmax;
     Istr.dstep(t).ymin = ymin;
     Istr.dstep(t).ymax = ymax;
    end
    Istr.z1 = zeros(m1*m2, 1);

    Ostr = most(Istr, mpopt);

    s = load(solnfile);

    t = 'dynamical system state (Z)';
    t_is(Ostr.results.Z, s.Z, 4, t);
else
    t_skip(1, 'requires CPLEX, Gurobi, MOSEK or quadprog');
end

% YorN = input('Play movie? (y/n) : ', 's');
% if strcmp(upper(YorN(1)), 'Y')
%     domovie;
% end

t_end

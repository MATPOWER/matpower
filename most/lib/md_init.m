function md = md_init
%MD_INIT  MOST Data structure Initialization
%   MD = MD_INIT
%
%   Creates an empty MOST Data struct (MD) with all fields required for
%   MOST, both input and output fields.

%   MOST
%   Copyright (c) 2010-2016, Power Systems Engineering Research Center (PSERC)
%   by Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Nacional de Colombia
%   and Ray Zimmerman, PSERC Cornell
%
%   This file is part of MOST.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/most for more info.


%%-----  start over  -----
md = struct('Delta_T', 1);
md.Storage.UnitIdx                    = [];
md.UC.CommitKey                       = [];   % If empty, no UC dimension; else, must contain key to be used
md.Storage.ExpectedTerminalStorageAim = [];   % expected terminal storage targe
md.Storage.ExpectedTerminalStorageMin = [];   % expected terminal storage LB
md.Storage.ExpectedTerminalStorageMax = [];   % expected terminal storage UB
md.Storage.rho                        = [];   % (i,t), varies bounding of storage dispatches between being
                                              % based on worst case (rho=1) vs. expected (rho=0) stored energy
                                              % at beginning of period
md.Storage.TerminalChargingPrice0     = [];   % applied to psc_tij0 (end-of-horizon terminal states)
md.Storage.TerminalDischargingPrice0  = [];   % applied to psd_tij0 (end-of-horizon terminal states)
md.Storage.TerminalChargingPriceK     = [];   % applied to psc_tijk (contingency terminal states)
md.Storage.TerminalDischargingPriceK  = [];   % applied to psd_tijk (contingency terminal states)


% %%-----  problem input data  -----
% md = struct('OpenEnded', 1);  % default = no terminal dispatch ramping constraints
% md.QCoordination  = 0;    % Create Qg variables for coordination if 1.
% md.TerminalPg     = [];   % dispatch to ramp to from the final period
% md.InitialPg      = [];   % dispatch to ramp from in the initial period
% 
% % used only when most.storage.terminal_target option is 1
% md.Storage.ExpectedTerminalStorageAim     = [];   % expected terminal storage targe
% md.Storage.ExpectedTerminalStorageMin     = [];   % expected terminal storage LB
% md.Storage.ExpectedTerminalStorageMax     = [];   % expected terminal storage UB
% 
% md.mpc.baseMVA    = [];
% md.mpc.bus        = [];
% md.mpc.gen        = [];
% md.mpc.branch     = [];
% md.mpc.gencost    = [];
% %md.offer(1).gencost                           = []; (deprecated)
% md.offer(1).PositiveActiveReservePrice        = [];   % (t)
% md.offer(1).PositiveActiveReserveQuantity     = [];
% md.offer(1).NegativeActiveReservePrice        = [];
% md.offer(1).NegativeActiveReserveQuantity     = [];
% md.offer(1).PositiveActiveDeltaPrice          = [];
% md.offer(1).NegativeActiveDeltaPrice          = [];
% %md.offer(1).PositiveReactiveReservePrice      = [];
% %md.offer(1).PositiveReactiveReserveQuantity   = [];
% %md.offer(1).NegativeReactiveReservePrice      = [];
% %md.offer(1).NegativeReactiveReserveQuantity   = [];
% %md.offer(1).PositiveReactiveDeltaPrice        = [];
% %md.offer(1).NegativeReactiveDeltaPrice        = [];
% md.offer(1).PositiveLoadFollowReservePrice    = [];
% md.offer(1).PositiveLoadFollowReserveQuantity = [];
% md.offer(1).NegativeLoadFollowReservePrice    = [];
% md.offer(1).NegativeLoadFollowReserveQuantity = [];
% md.RampWearCostCoeff = [];    % (i,t) note different scheme!
%                                 % the first column is the cost from initial
%                                 % state to t=1; if there's a terminal
%                                 % state, then there must be nt+1 columns
% md.UC.CommitSched                     = [];   % if UC: solution on output; if not UC (mdi.CommitKey empty) then
%                                                 % must contain UC status
%                                                 % (i, t)
% md.UC.CommitKey                       = [];   % If empty, no UC dimension; else, must contain key to be used
%                                                 % for UC problem as follows:
%                                                 % Must run        : 2
%                                                 % Available for UC: 0, 1
%                                                 % Offline:        : -1
%                                                 % (i, t)
% md.UC.InitialState                    = [];   % If positive, number of uptime periods;
%                                                 % if negative, number of downtime periods at t = 0
%                                                 % (ng)
% md.UC.MinUp                           = [];   % Minimum uptime (i)
% md.UC.MinDown                         = [];   % Minimum downtime (i)
% md.UC.CyclicCommitment                = [];   % 1 if Commit restrictions roll over
% md.UC.c00                             = [];   % (i,t)
% md.Delta_T                            = 1;    % length of each period in hours
% md.Storage.UnitIdx                    = [];
% md.Storage.MinStorageLevel            = [];   % (i,t) These are applied to the
% md.Storage.MaxStorageLevel            = [];   % (i,t) bounds.
% md.Storage.InitialStorage             = [];
% md.Storage.InitialStorageLowerBound   = [];   % initial s- value, or (if ForceCyclicStorage) lower bound on s0
% md.Storage.InitialStorageUpperBound   = [];   % initial s+ value, or (if ForceCyclicStorage) upper bound on s0
% md.Storage.OutEff                     = [];   % (i,t), output efficency (if empty defaults to 1)
% md.Storage.InEff                      = [];   % (i,t), output efficency (if empty defaults to 1)
% md.Storage.InitialStorageCost         = [];   % this is the cost of s0 (not adjusted for output efficiency)
% md.Storage.TerminalStoragePrice       = [];   % applied to psc_tij0, psd_tij0 (non-terminal states)
% md.Storage.TerminalChargingPrice0     = [];   % applied to psc_tij0 (end-of-horizon terminal states)
% md.Storage.TerminalDischargingPrice0  = [];   % applied to psd_tij0 (end-of-horizon terminal states)
% md.Storage.TerminalChargingPriceK     = [];   % applied to psc_tijk (contingency terminal states)
% md.Storage.TerminalDischargingPriceK  = [];   % applied to psd_tijk (contingency terminal states)
% md.Storage.LossFactor                 = [];   % (i,t), fraction of storage lost per hour
%                                                 % if empty, defaults to lossless (all zeros)
% md.Storage.rho                        = [];   % (i,t), varies bounding of storage dispatches between being
%                                                 % based on worst case (rho=1) vs. expected (rho=0) stored energy
%                                                 % at beginning of period
% % Note: MinStorageLevel, MaxStorageLevel, InEff, OutEff, LossFactor and rho
% %       are optionally expanded automatically in most() from scalar,
% %       ns x 1, or 1 x nt to ns x nt matrix.
% 
% md.idx.nt                             = [];
% md.cont(1,1).contab                   = [];
% md.tstep(1).OpCondSched(1).tab        = [];
% md.tstep(1).TransMat = [];    % tstep(t).TransMat(j(t), j(t-1)).
%                                 % Note that for cyclic or terminal state
%                                 % problems, data for t=nt+1 is needed for
%                                 % this matrix!
% md.tstep(1).TransMask = [];   % Same format as TransMat, mask indicating
%                                 % whether to include transition in ramp reserve.
% md.tstep(1).Li = sparse(0,0); % These 8 for computing expected storage states
% md.tstep(1).Lf = sparse(0,0);
% md.tstep(1).Mg = sparse(0,0);
% md.tstep(1).Mh = sparse(0,0);
% md.tstep(1).Ng = sparse(0,0);
% md.tstep(1).Nh = sparse(0,0);
% md.tstep(1).G  = sparse(0,0);
% md.tstep(1).H  = sparse(0,0);
% md.tstep(1).E  = sparse(0,0); % To compute expected injections in t-th period
% 
% md.dstep(1).A = sparse(0,0);
% md.dstep(1).B = sparse(0,0);
% md.dstep(1).C = sparse(0,0);
% md.dstep(1).D = sparse(0,0);
% md.dstep(1).zmax = [];
% md.dstep(1).zmin = [];
% md.dstep(1).ymax = [];
% md.dstep(1).ymin = [];
% md.idx.ntds = []; % Number of periods in the dynamical system horizon
% md.z1 = [];       % Initial state (t=1)
% 
% %%-----  internally created data:  -----
% %% (1) Indexing mechanism
% md.idx.nj         = [];
% md.idx.nc         = [];
% md.idx.nb         = [];
% md.idx.nb_total   = [];
% md.idx.ng         = [];
% md.idx.nf_total   = [];
% md.idx.ns         = [];
% md.idx.ns_total   = [];
% md.idx.nzds       = [];       % size of state vector for dynamical system
% md.idx.nyds       = [];
% md.idx.ntramp     = [];       % number of periods of load following reserves
% % md.idx.thbas      = [];
% % md.idx.thend      = [];
% % md.idx.pbas       = [];
% % md.idx.pend       = [];       % (t,j,k)
% % md.idx.dppbas     = [];       % (t,j,k)
% % md.idx.dppend     = [];       % (t,j,k)
% % md.idx.dpmbas     = [];       % (t,j,k)
% % md.idx.dpmend     = [];       % (t,j,k)
% % md.idx.ybas       = [];       % (t,j,k)
% % md.idx.yend       = [];       % (t,j,k)
% % md.idx.pcbas      = [];       % (t)
% % md.idx.pcend      = [];       % (t)
% % md.idx.rppbas     = [];       % (t)
% % md.idx.rppend     = [];       % (t)
% % md.idx.rpmbas     = [];       % (t)
% % md.idx.rpmend     = [];       % (t)
% % md.idx.pscbas     = [];       % (t,j,k)
% % md.idx.pscend     = [];       % (t,j,k)
% % md.idx.psdbas     = [];       % (t,j,k)
% % md.idx.psdend     = [];       % (t,j,k)
% % md.idx.rrpbas     = [];       % (t)
% % md.idx.rrpend     = [];       % (t)
% % md.idx.rrmbas     = [];       % (t)
% % md.idx.rrmend     = [];       % (t)
% % md.idx.spbas      = [];       % (t)
% % md.idx.spend      = [];       % (t)
% % md.idx.smbas      = [];       % (t)
% % md.idx.smend      = [];       % (t)
% % md.idx.s0bas      = [];       % (1:ns)
% % md.idx.s0end      = [];       % (i:ns)
% % md.idx.ubas       = [];       % (t)
% % md.idx.uend       = [];       % (t)
% % md.idx.wbas       = [];       % (t)
% % md.idx.wend       = [];       % (t)
% % md.idx.vbas       = [];       % (t)
% % md.idx.vend       = [];       % (t)
% % md.idx.qbas       = [];       % (t)
% % md.idx.qend       = [];       % (t)
% % md.idx.netbas     = [];       % (t,j,k)
% % md.idx.netend     = [];       % (t,j,k)
% % md.idx.lfbas      = [];       % (t,j,k)
% % md.idx.lfend      = [];       % (t,j,k)
% % md.idx.lstibas    = [];       % (t,j,k)
% % md.idx.lstiend    = [];       % (t,j,k)
% % md.idx.Aybas      = [];       % (t,j)
% % md.idx.Ayend      = [];       % (t,j)
% % md.idx.lc1bas     = [];       % (t)
% % md.idx.lc1end     = [];       % (t)
% % md.idx.lc2bas     = [];       % (t)
% % md.idx.lc2end     = [];       % (t)
% % md.idx.lc3bas     = [];       % (t,j,k)
% % md.idx.lc3end     = [];       % (t,j,k)
% % md.idx.lc5bas     = [];       % (t,j,k)
% % md.idx.lc5end     = [];       % (t,j,k)
% % md.idx.lc6bas     = [];       % (t,j,k)
% % md.idx.lc6end     = [];       % (t,j,k)
% % md.idx.lc7bas     = [];       % (t,j,k)
% % md.idx.lc7end     = [];       % (t,j,k)
% % md.idx.lc8bas     = [];       % (t,j,k)
% % md.idx.lc8end     = [];       % (t,j,k)
% % md.idx.lc9bas     = [];       % (t,j,k)
% % md.idx.lc9end     = [];       % (t,j,k)
% % md.idx.lc10bas    = [];       % (t,j,k)
% % md.idx.lc10end    = [];       % (t,j,k)
% % md.idx.lc21bas    = [];       % (t,jt,j(t+1))
% % md.idx.lc21end    = [];       % (t,jt,j(t+1))
% % md.idx.lc22bas    = [];       % (t)
% % md.idx.lc22end    = [];       % (t)
% % md.idx.lc23bas    = [];       % (t,jt,j(t+1))
% % md.idx.lc23end    = [];       % (t,jt,j(t+1))
% % md.idx.lc24bas    = [];       % (t)
% % md.idx.lc24end    = [];       % (t)
% % md.idx.lc31bas    = [];       % (t)
% % md.idx.lc31end    = [];       % (t)
% % md.idx.lc32bas    = [];       % (t)
% % md.idx.lc32end    = [];       % (t)
% % md.idx.lc33bas    = [];       % (t)
% % md.idx.lc33end    = [];       % (t)
% % md.idx.lc34bas    = [];       % (t,j)
% % md.idx.lc34end    = [];       % (t,j)
% % md.idx.lc35bas    = [];       % (t,j)
% % md.idx.lc35end    = [];       % (t,j)
% % md.idx.lc36bas    = [];       % (t,j)
% % md.idx.lc36end    = [];       % (t,j)
% % md.idx.lc37bas    = [];       % (t,j)
% % md.idx.lc37end    = [];       % (t,j)
% % md.idx.lc38bas    = [];       %
% % md.idx.lc38end    = [];       %
% % md.idx.lc40bas    = [];       % (t)
% % md.idx.lc40end    = [];       % (t)
% % md.idx.lc41bas    = [];       % (t)
% % md.idx.lc41end    = [];       % (t)
% % md.idx.lc50bas    = [];       % (t)
% % md.idx.lc50end    = [];       % (t)
% % md.idx.lc51bas    = [];       % (i,t)
% % md.idx.lc51end    = [];       % (i,t)
% % md.idx.lc52bas    = [];       % (i,t)
% % md.idx.lc52end    = [];       % (i,t)
% % md.idx.lc53bas    = [];       % (t,j,k)
% % md.idx.lc53end    = [];       % (t,j,k)
% % md.idx.lc54bas    = [];       % (t,j,k)
% % md.idx.lc54end    = [];       % (t,j,k)
% % md.idx.lc55bas    = [];       % (t,j,k)
% % md.idx.lc55end    = [];       % (t,j,k)
% % md.idx.lc56bas    = [];       % (t,j,k)
% % md.idx.lc56end    = [];       % (t,j,k)
% 
% md.flow(1,1,1).mpc = md.mpc;
% md.DCMODEL              = [];   % DC flow used to model the network as opposed
%                                 % to simple generation = demand constraint
%                                 % (set via mpopt.most.dc_line)
% md.SecurityConstrained  = [];   % contingencies cases considered
%                                 % (set via mpopt.most.security_constraints and
%                                 %  presence (or not) of contingency data)
% md.Storage.ForceExpectedTerminalStorage   = [];     % flag, 0 or 1, terminal storage target included
% md.Storage.ForceCyclicStorage             = [];     % 1 = includes cyclic constraint (initial storage
%                                                     % is a var = final expected storage), 0 = does not include
% md.UC.run               = [];   % 1 = run the unit commitment, 0 = don't
% md.alpha                = [];   % defines when during period contingencies happen
% 
% %% (2) QP problem
% md.QP.A       = sparse(0,0);
% md.QP.l       = [];
% md.QP.u       = [];
% md.QP.xmin    = [];
% md.QP.xmax    = [];
% md.QP.vtype   = '';
% md.QP.H       = sparse(0,0);
% md.QP.C1      = [];
% md.QP.c1      = [];
% md.QP.C       = [];
% md.QP.Cfstor  = [];
% md.CoordCost.Huser    = sparse(0,0);
% md.CoordCost.Cuser    = [];
% md.CoordCost.cuser    = [];
% md.QP.x           = [];
% md.QP.f           = [];
% md.QP.exitflag    = [];
% md.QP.output      = [];
% md.QP.lambda      = [];
% md.QP.opt         = [];
% 
% %%-----  result data  -----
% md.results.f      = [];
% md.results.Pc     = [];
% md.results.Rpp    = [];
% md.results.Rpm    = [];
% md.results.Rrp    = [];
% md.results.Rrm    = [];
% md.results.Sp     = [];
% md.results.Sm     = [];
% md.results.GenPrices          = [];
% md.results.CondGenPrices      = [];
% md.results.RrpPrices          = [];
% md.results.RrmPrices          = [];
% md.results.ExpectedRampCost   = [];
% md.results.SetupTime          = [];
% md.results.SolveTime          = [];
% md.results.Z                  = [];
% md.results.Y                  = [];
% 
% md.CostWeights    = [];   % (k,j,t) !!!! NOTE order! So that (:,:,t)
%                             % refers to t-th period
% md.CostWeightsAdj = [];   % (k,j,t) !!!! NOTE order! So that (:,:,t)
%                             % refers to t-th period
% md.StepProb       = [];   % (t)  - probability of making it to the t-th step
% md.Storage.ExpectedStorageState       = [];
% md.Storage.ExpectedStorageDispatch    = [];

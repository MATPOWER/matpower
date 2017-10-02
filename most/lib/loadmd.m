function md = loadmd(mpci, transmati, xgdi, storagei, contabi, profilesi, trajdatai)
%LOADMD   Loads all required data and constructs MD for MOST.
%
%   MD = LOADMD(MPC)
%   MD = LOADMD(MPC, TRANSMAT)
%   MD = LOADMD(MPC, TRANSMAT, XGD)
%   MD = LOADMD(MPC, TRANSMAT, XGD, SD)
%   MD = LOADMD(MPC, TRANSMAT, XGD, SD, CONTAB)
%   MD = LOADMD(MPC, TRANSMAT, XGD, SD, CONTAB, PROFILES)
%   MD = LOADMD(MPC, TRANSMAT, XGD, SD, CONTAB, PROFILES, TRAJDATA)
%   MD = LOADMD(MPC, NT, ...)
%
%   All inputs can have the format described below or can be given as
%   strings containing the name of the MAT-file or M-file that contains the
%   data in the proper format.
%
%   Inputs:
%       MPC:       a standard MATPOWER case struct, optionally with
%                  additional fields such as 'genfuel' and 'i<type>'
%                  NOTE: Bus numbers must be consecutive beginning at 1
%                        (i.e. internal ordering).
%       TRANSMAT:  (optional) NT dimensional cell array of matrices, where
%                  TRANSMAT{t} is an NJ(t) x NJ(t-1) matrix containing the
%                  transition probabilities from period t-1 to period t. The
%                  first element TRANSMAT{1} is a column vector of transition
%                  probabilities from period 0 (NJ(0) = 1) to period 1. For
%                  deterministic cases, TRANSMAT can be specified simply as
%                  an integer NT (number of periods), which gets expanded
%                  internally to a cell array of 1's. Default value is 1.
%       XGD:       (optional) xGenData struct, see LOADXGENDATA for details.
%       SD:        (optional) StorageData struct, see LOADSTORAGEDATA for
%                  details.
%       CONTAB:    (optional) contingency table with master set of
%                  contingencies used for security throughout entire horizon
%       PROFILES:  (optional) a struct array of Profiles (see IDX_PROFILE
%                  and APPLY_PROFLE for details), specifying changes in the
%                  system and operational conditions, xGenData, StorageData
%                  and/or contingencies across time periods and scenarios.
%                  Alternatively, for backward compatibility with the
%                  older centroids format, PROFILES can take the form of
%                  a struct with the following fields:
%                     .wind:  3-dim array (NT x NJ_MAX x num wind sites)
%                     .load:  3-dim array (NT x NJ_MAX x num load zones)
%                  If empty, no changes are made across time.
%       TRAJDATA:  (optional) struct array in same format as PROFILES
%                  specifying a set of trajectories for a stage 2 solution.
%                  The J dimension of the profile in this case indexes the
%                  trajectories.

% Created by Daniel Munoz-Alvarez (2/28/2013)

%   MOST
%   Copyright (c) 2013-2016, Power Systems Engineering Research Center (PSERC)
%   by Daniel Munoz-Alvarez and Ray Zimmerman, PSERC Cornell
%
%   This file is part of MOST.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/most for more info.

if nargin < 7
    trajdatai = [];
    if nargin < 6
        profilesi = [];
        if nargin < 5
            contabi = [];
            if nargin < 4
                storagei = [];
                if nargin < 3
                    xgdi = [];
                    if nargin < 2
                        transmati = [];
                        if nargin < 1
                            error('loadmd: MPC is a mandatory argument');
                        end
                    end
                end
            end
        end
    end
end

define_constants;
[CT_LABEL, CT_PROB, CT_TABLE, CT_TBUS, CT_TGEN, CT_TBRCH, ...
    CT_TAREABUS, CT_TAREAGEN, CT_TAREABRCH, CT_ROW, CT_COL, CT_CHGTYPE, ...
    CT_REP, CT_REL, CT_ADD, CT_NEWVAL, CT_TLOAD, CT_TAREALOAD, ...
    CT_LOAD_ALL_PQ, CT_LOAD_FIX_PQ, CT_LOAD_DIS_PQ, CT_LOAD_ALL_P, ...
    CT_LOAD_FIX_P, CT_LOAD_DIS_P, CT_TGENCOST, CT_TAREAGENCOST, ...
    CT_MODCOST_F, CT_MODCOST_X] = idx_ct;
[PR_REP, PR_REL, PR_ADD, PR_TCONT, PR_TYPES, PR_TMPCD,...
   PR_TXGD, PR_TCTD, PR_TSTGD, PR_CHGTYPES] = idx_profile;

% (A) load mandatory arg MPC

% (A.1) load MATPOWER case
mpc = loadcase(mpci);


% (B) read other arguments into struct

% (B.1) transmat
if isempty(transmati)
    transmat = {1};         % default is single period deterministic
elseif isscalar(transmati) && isnumeric(transmati)
    transmat = transmati;   % specifies NT
else
    type = {'cell'};
    transmat = loadgenericdata(transmati, type);
    if isempty(transmat)
        error('loadmd: must provide a non-empty TRANSMAT cell array')
    end
end

% (B.2) xGenData
type = {'struct'};
fields = {...
          'CommitSched', ...
          'InitialPg', ...
          'RampWearCostCoeff', ...
          'PositiveActiveReservePrice', ...
          'PositiveActiveReserveQuantity', ...
          'NegativeActiveReservePrice', ...
          'NegativeActiveReserveQuantity', ...
          'PositiveActiveDeltaPrice', ...
          'NegativeActiveDeltaPrice', ...
          'PositiveLoadFollowReservePrice', ...
          'PositiveLoadFollowReserveQuantity', ...
          'NegativeLoadFollowReservePrice', ...
          'NegativeLoadFollowReserveQuantity', ...
          };
xgd = loadgenericdata(xgdi, type, fields);


% (B.3) storage
type = {'struct'};
fields = {...
          'UnitIdx';...
%           'ExpectedTerminalStorageAim';... optional
%           'ExpectedTerminalStorageMin';... optional
%           'ExpectedTerminalStorageMax';... optional
          'InitialStorage';...
%           'InitialStorageLowerBound';... optional
%           'InitialStorageUpperBound';... optional
%           'InitialStorageCost';... optional
          'TerminalStoragePrice';...
%           'TerminalChargingPrice0';... optional
%           'TerminalDischargingPrice0';... optional
%           'TerminalChargingPriceK';... optional
%           'TerminalDischargingPriceK';... optional
%           'ExpectedStorageState';... optional
%           'ExpectedStorageDispatch';... optional
          'MinStorageLevel';...
          'MaxStorageLevel';...
%           'OutEff';... optional
%           'InEff';... optional
%           'LossFactor';... optional
%           'rho';... optional
          };
storage = loadgenericdata(storagei, type, fields);


% (B.4) contab
type = {'array'};
contab = loadgenericdata(contabi, type);


% (B.5) profiles
type = {'struct'};
fields = {'type','table','rows','col','chgtype','values'};

% Load profiles. If fails to find profiles' fields, then checks whether the
% provided input is in a centroids format. This is obsolete though since
% centroids only supports changes across time for load and wind. Keeping
% just for backard compatibility.
try
    profiles = loadgenericdata(profilesi, type, fields);
catch exception
    if strcmp(exception.identifier,'loadgenericdata:missingfield')
        fields = {'wind', 'load'};
        centroids = loadgenericdata(profilesi, type, fields);
        if size(centroids.load,3) > 1 % if 3rd-dim = several zones
            table = PR_TAREALOAD;
            rows = (1:size(centroids.load,3))';
        else
            table = PR_TLOAD;
            rows = 0;
        end
        profiles = [ centroids2profile(centroids.wind, PR_TGEN, mpc.iwind, PMAX, PR_REL) ;...
                     centroids2profile(centroids.load, table, rows, CT_LOAD_ALL_P, PR_REL) ;
                   ];
    else
        rethrow(exception);
    end
end


% (B.6) trajdata
type = {'struct'};
fields = {'type','table','rows','col','chgtype','values'};

% Load trajdata works as load profile. See details above regarding the
% backward compatibility to use centroids format.
try
    trajdata = loadgenericdata(trajdatai, type, fields);
catch exception
    if strcmp(exception.identifier,'loadgenericdata:missingfield')
        fields = {'wind', 'load'};
        trajdatacentroids = loadgenericdata(trajdatai, type, fields);
        if size(trajdatacentroids.load,3) > 1 % if 3rd-dim = several zones
            table = PR_TAREALOAD;
            rows = (1:size(trajdatacentroids.load,3))';
        else
            table = PR_TLOAD;
            rows = 0;
        end
        trajdata = [ centroids2profile(trajdatacentroids.wind, PR_TGEN, mpc.iwind, PMAX, PR_REL) ;...
                     centroids2profile(trajdatacentroids.load, table, rows, CT_LOAD_ALL_P, PR_REL) ;
                   ];
    else
        rethrow(exception);
    end
end
if ~isempty(trajdata) % trajdatai not provided causes info regarding second stg not to be loaded
    second = true;
else
    second = false;
end


% (C) verify consistency of inputs and truncate them to run algorithm when possible

% (C.0) Initialize md struct
md = md_init;


% (C.1) Obtain basic dimensions from the most restrictive input data
if iscell(transmat)
    nt = length(transmat);
else                            % expand NT to full deterministic TRANSMAT
    nt = transmat;
    transmat = ones(1, nt);
    transmat = mat2cell(transmat, 1, transmat);
end

nj = zeros(nt,1);
for t = 1:nt
    nj(t) = size(transmat{t},1);
end
nj_max = max(nj);

nb = size(mpc.bus,1);           % num buses
ng = size(mpc.gen,1);           % num gen
if isfield(mpc, 'iwind')
    nw = size(mpc.iwind,1);     % num wind gen
else
    nw = 0;                     % num wind gen
end
if isfield(mpc, 'iess')
    ns = size(mpc.iess,1);      % num ess gen
else
    ns = 0;                     % num ess gen
end
nprof = size(profiles,1);       % num of profiles (wind profile, load profile, ice storage profile, etc.)

if second
    nprof2 = size(trajdata,1);              % Determined as num of profile elements (structs) there are in trajdata
    ntraj = size(trajdata(1).values,2);     % Determined as num of trajectories that the first trajectory profile has
end


% (C.2) Store basic dimensions not imposed implicitly by inputs
    %     themselves in most(.)
md.idx.nt = nt;


% (C.3) Verify that obtained basic dimensions are valid
if any(nj < 1)  % verify nj dimensions
    error('loadmd: number of scenarios nj must be at least 1');
end

if nt < 1       % need at least one time period
    error('loadmd: number of time periods nt must be at least 1');
end


% (C.4) Look for empty input data and assign well-formatted default data
    % which should be assigned in case the corresponding variable is empty.

% (C.4.0) mpc
%% check that bus numbers are equal to indices to bus (one set of bus numbers)
if any(mpc.bus(:, BUS_I) ~= (1:nb)')
    error('loadmd: buses must be numbered consecutively in MPC.bus matrix; use ext2int() to convert to internal ordering')
end

% (C.4.1) transmat
    % handled previously so dimensions could be set in (C.1)

% (C.4.2) xGenData
if isempty(xgd)     % assign the default input
    xgd = loadxgendata([], mpc);
end

% (C.4.3) storage
if ns == 0  % No ess units
    storage = [];   % storage struct ignored if provided
else
    if isempty(storage)
        error('loadmd: storage struct cannot be empty when MPC contains storage units');
    end
    % let most() assign defaults
    if ~isfield(storage,'OutEff')
        storage.OutEff = [];
    end
    if ~isfield(storage,'InEff')
        storage.InEff = [];
    end
    if ~isfield(storage,'LossFactor')
        storage.LossFactor = [];
    end
    if ~isfield(storage,'rho')
        storage.rho = [];
    end
end

% (C.4.4) contab

% (C.4.5) profiles
    % Optional: An empty profiles struct array means no differences in the
    % operating conditions across time periods or scenarios.

% (C.4.6) options

% (C.4.7) trajdata
    % Optional, if not given (empty) then it is not used


% (C.5) Verify consistency of input dimensions

% (C.5.1) transmat
if length(transmat) < nt   % verify transmat dim is nt (redundant since nt is defined with transmat)
    error('loadmd: nt (= %d) has been modified and no longer matches transmat time dimension (%d)', nt, length(transmat));
else        % verify consistency of transmat with nj(t) (redundant since nj is defined with transmat)
    if size(transmat{1}, 1) ~= nj(1)
        error('loadmd: # of scenarios in t=1, nj(1) = %d, does not match # of rows in transmat{1} (%d)', nj(1), size(transmat{1}, 1));
    end
    if nt > 1
        for t = 2:nt
            if any(size(transmat{t}) ~= [nj(t) nj(t-1)])
                error('loadmd: dimensions of transmat{%d} (%d x %d) inconsistent with nj(%d) x nj(%d) (%d x %d)', ...
                    t, size(transmat{t}, 1), size(transmat{t}, 2), t, t-1, nj(t), nj(t-1));
            end
        end
    end
end

% (C.5.2) xGenData
if any( [size(xgd.CommitSched, 1) ...
         size(xgd.InitialPg, 1) ...
         size(xgd.RampWearCostCoeff, 1) ...
         size(xgd.PositiveActiveReserveQuantity, 1) ...
         size(xgd.NegativeActiveReservePrice, 1) ...
         size(xgd.NegativeActiveReserveQuantity, 1) ...
         size(xgd.PositiveActiveDeltaPrice, 1) ...
         size(xgd.NegativeActiveDeltaPrice, 1) ...
         size(xgd.PositiveLoadFollowReservePrice, 1) ...
         size(xgd.PositiveLoadFollowReserveQuantity, 1) ...
         size(xgd.NegativeLoadFollowReservePrice, 1) ...
         size(xgd.NegativeLoadFollowReserveQuantity, 1) ...
         ] ~= ng )
     error('loadmd: 1st dimension inconsistency in field of xGenData')
end

% (C.5.3) storage
if ~isempty(storage)
    if length(storage.UnitIdx) ~= length(mpc.iess)
        error('loadmd: dimensions of storage.UnitIdx (%d) and mpc.iess (%d) must match', length(storage.UnitIdx), length(mpc.iess));
    end
end

% (C.5.4) contab
size_contab = size(contab);
if ~isempty(contab) && (length(size_contab) ~= 2 || size_contab(2) ~= 7)
    error('loadmd: contab must be matrix with 7 columns');
end

% (C.5.5) profiles
for p = 1:nprof
    if ~isnumeric(profiles(p).values)
        error('loadmd: profiles(%d).values is required to be a numeric array', p);
    end
    if length(size(profiles(p).values)) > 3
        error('loadmd: profiles(%d).values must have no more than 3 dimensions', p);
    end
    if ~( size(profiles(p).values, 1) == 1 ||...
          size(profiles(p).values, 1) >= nt )
        error('loadmd: time dimension of profiles(%d).values (%d) must be 1 or >= nt = %d', p, size(profiles(p).values, 1), nt);
    end
    if ~( size(profiles(p).values, 2) == 1 ||...
          size(profiles(p).values, 2) == nj_max )
        error('loadmd: scenarios dimension of profiles(%d).values (%d) must be 1 or nj_max = %d', p, size(profiles(p).values, 2), nj_max);
    end
    if ~( size(profiles(p).values, 3) == 1 ||...
          size(profiles(p).values, 3) == length(profiles(p).rows) )
        error('loadmd: 3rd dimension of profiles(%d).values (%d) must be 1 or %d (length of profiles(%d).rows)', p, size(profiles(p).values, 3), length(profiles(p).rows), p);
    end
    if ~any(strcmp(profiles(p).type, PR_TYPES))
        error('loadmd: profiles(%d).type not supported', p);
    end
    if ~(isscalar(profiles(p).table) || ischar(profiles(p).table))
        error('loadmd: profiles(%d).table must be either scalar or string', p);
    end
    if ~isvector(profiles(p).rows)
        error('loadmd: profiles(%d).rows must be a vector', p);
    end
    if ~isscalar(profiles(p).col)
        error('loadmd: profiles(%d).col must be scalar', p);
    end
    if ~any(profiles(p).chgtype == PR_CHGTYPES)
        error('loadmd: profiles(%d).chgtype not supported', p);
    end
end

% (C.5.6) options

% (C.5.7) trajdata
if second
    for p = 1:nprof2
        if ~isnumeric(trajdata(p).values)
            error('loadmd: trajdata(%d).values is required to be a numeric array', p);
        end
        if length(size(trajdata(p).values)) > 3
            error('loadmd: trajdata(%d).values must have no more than 3 dimensions', p);
        end
        if ~( size(trajdata(p).values, 1) == 1 ||...
              size(trajdata(p).values, 1) >= nt )
            error('loadmd: time dimension of trajdata(%d).values (%d) must be 1 or >= nt = %d', p, size(trajdata(p).values, 1), nt);
        end
        if ~( size(trajdata(p).values, 2) == 1 ||...
              size(trajdata(p).values, 2) == ntraj )
            error('loadmd: scenarios dimension of trajdata(%d).values (%d) must be 1 or ntraj = %d', p, size(trajdata(p).values, 2), ntraj);
        end
        if ~( size(trajdata(p).values, 3) == 1 ||...
              size(trajdata(p).values, 3) == length(trajdata(p).rows) )
            error('loadmd: 3rd dimension of trajdata(%d).values (%d) must be 1 or %d (length of trajdata(%d).rows)', p, size(trajdata(p).values, 3), length(trajdata(p).rows), p);
        end
        if ~any(strcmp(trajdata(p).type, PR_TYPES))
            error('loadmd: trajdata(%d).type not supported', p);
        end
        if ~(isscalar(trajdata(p).table) || ischar(trajdata(p).table))
            error('loadmd: trajdata(%d).table must be either scalar or string', p);
        end
        if ~isvector(trajdata(p).rows)
            error('loadmd: trajdata(%d).rows must be a vector', p);
        end
        if ~isscalar(trajdata(p).col)
            error('loadmd: trajdata(%d).col must be scalar', p);
        end
        if ~any(trajdata(p).chgtype == PR_CHGTYPES)
            error('loadmd: trajdata(%d).chgtype not supported', p);
        end
    end
end


% (C.6) Truncate all input data to make consistent basic dimensions
    % Not necessary for transmat (since it defines nt)

% (C.6.1) transmat
transmat = transmat(1:nt);  % redundant, since nt is defined by transmat
% (C.6.2) xGenData
% (C.6.3) storage
% (C.5.4) contab
% (C.6.5) profiles
for p = 1:nprof
    profiles(p).values = profiles(p).values(1:nt,:,:);
end
% (C.6.6) options
% (C.6.7) trajdata
if second
    for p = 1:nprof2
        trajdata(p).values = trajdata(p).values(1:nt,1:ntraj,:);
    end
end

% (C.7) Put all necessary info into md struct in order to run most()
%       Here, all the profiles are "expanded" and applied to data before
%       the latter is actually assigned to md.
%       Moreover, all info required to run mpsopfl2 is loaded into md if
%       second == 1.

% (C.7.0) transmat: assign transition probability matrices to md struct
for t = 1:nt
    md.tstep(t).TransMat = transmat{t};
end

% (C.7.1) mpc: assign mpc to md struct
md.mpc = mpc;

% (C.7.2) profiles: build operation conditions contab per period per scenario

% (C.7.2.1) transform profiles of type mpcData into contabs
optab = cell(nt,nj_max);
for p = 1:nprof
    if strcmp(profiles(p).type, 'mpcData')
%       Profiles of type mpcData need to be passed along to most()
%       as contingency-like tables describing the operation
%       conditions for each time period and each scenario. Those
%       contabs are stored in the OpCondSched struct array.
        optab = apply_profile(profiles(p), optab);
    end
end

% (C.7.2.2) store 'operation conditions' contabs into md
for t = 1:nt
    for j = 1:nj(t)
        md.tstep(t).OpCondSched(j).tab = optab{t,j};
    end
end
    
% (C.7.2.3) store profiles into md perhaps for results
% md.profiles = profiles;   % For ploting purposes only (need wind profile in mpsopf_plots1)

% (C.7.3) storage: superimpose given storage struct into initialized
%                  storage struct

% (C.7.3.1) assign fields given in 'storage' input to initialized 'Storage' struct
if ~isempty(storage)
    fields = { ...
        'UnitIdx', ...
        'ExpectedTerminalStorageAim', ...
        'ExpectedTerminalStorageMin', ...
        'ExpectedTerminalStorageMax', ...
        'InitialStorage', ...
        'InitialStorageLowerBound', ...
        'InitialStorageUpperBound', ...
        'InitialStorageCost', ...
        'TerminalStoragePrice', ...
        'TerminalChargingPrice0', ...
        'TerminalDischargingPrice0', ...
        'TerminalChargingPriceK', ...
        'TerminalDischargingPriceK', ...
        'ExpectedStorageState', ...
        'ExpectedStorageDispatch', ...
        'MinStorageLevel', ...
        'MaxStorageLevel', ...
        'OutEff', ...
        'InEff', ...
        'LossFactor', ...
        'rho', ...
    };
    for f = 1:length(fields)
        ff = fields{f};
        if isfield(storage, ff)
            md.Storage.(ff) = storage.(ff);
        end
    end
end

% (C.7.3.2) apply storage profiles
if ns > 0
    for p = 1:nprof
        if strcmp(profiles(p).type, 'StorageData')
            md.Storage = apply_profile(profiles(p), md.Storage, ns);
        end
    end
end

% (C.7.4) options: assign options related variables

% (C.7.5) xGenData: assign xgd related variables
for p = 1:nprof
    if strcmp(profiles(p).type, 'xGenData')
        xgd = apply_profile(profiles(p), xgd, ng);
    end
end

if isfield(xgd, 'CommitKey') && ~isempty(xgd.CommitKey)
    UC = 1;
else
    UC = 0;
end

%% fields with no time dimension
%% in md
md.InitialPg = xgd.InitialPg;
if isfield(xgd, 'TerminalPg')       %% optional, deprecated
    md.TerminalPg = xgd.TerminalPg;
end
%% in md.UC
if UC
    md.UC.InitialState    = xgd.InitialState;
    md.UC.MinUp           = xgd.MinUp;
    md.UC.MinDown         = xgd.MinDown;
end

%% fields in md
fields = {'RampWearCostCoeff'};
for f = 1:length(fields)
    ff = fields{f};
    if size(xgd.(ff), 2) == 1       %% xGenData contains 1 column
        for t = 1:nt
            md.(ff)(:,t)   = xgd.(ff);
        end
    elseif size(xgd.(ff), 2) == nt  %% xGenData contains nt columns
        md.(ff)            = xgd.(ff);
    else
        error('loadmd: number of columns in XGD.%s must be 1 or %d', ff, nt);
    end
end

%% fields in md.UC
fields = {'CommitSched'};
if UC                               %% optional
    fields = {fields{:}, 'CommitKey'};
end
for f = 1:length(fields)
    ff = fields{f};
    if size(xgd.(ff), 2) == 1       %% xGenData contains 1 column
        for t = 1:nt
            md.UC.(ff)(:,t)   = xgd.(ff);
        end
    elseif size(xgd.(ff), 2) == nt  %% xGenData contains nt columns
        md.UC.(ff)            = xgd.(ff);
    else
        error('loadmd: number of columns in XGD.UC.%s must be 1 or %d', ff, nt);
    end
end

%% fields in offer
fields = {
    'PositiveActiveReservePrice', 'PositiveActiveReserveQuantity', ...
    'NegativeActiveReservePrice', 'NegativeActiveReserveQuantity', ...
    'PositiveActiveDeltaPrice', 'NegativeActiveDeltaPrice', ...
    'PositiveLoadFollowReservePrice', 'PositiveLoadFollowReserveQuantity', ...
    'NegativeLoadFollowReservePrice', 'NegativeLoadFollowReserveQuantity', ...
};
for f = 1:length(fields)
    ff = fields{f};
    if size(xgd.(ff), 2) == 1       %% xGenData contains 1 column
        for t = 1:nt
            md.offer(t).(ff)  = xgd.(ff);
        end
    elseif size(xgd.(ff), 2) == nt  %% xGenData contains nt columns
        for t = 1:nt
            md.offer(t).(ff)  = xgd.(ff)(:, t);
        end
    else
        error('loadmd: number of columns in XGD.%s must be 1 or %d', ff, nt);
    end
end

% (C.7.6) contab: assign contingency tables per period per scenario to md struct
ct_subset = ones(nt, nj_max, size(contab,1)); % subset of contingencies
%                                               in contab to be apply 
%                                               per scenario per time
%                                               period.
for p = 1:nprof
    if strcmp(profiles(p).type, 'ContingencyData')
        ct_subset = apply_profile(profiles(p), ct_subset, ncont);
    end
end

% Assign subset of contingencies indicated by profile ct_subset
for t = 1:nt
    for j = 1:nj(t)
        md.cont(t,j).contab = contab(squeeze(ct_subset(t,j,:)) == 1,:);
    end
end


% (C.7.7) trajdata: second stage related data
%       Although the data necessary for running the second stg is
%       loaded here, a first stage still needs to be ran before the
%       second stg can actually be ran since results from the former
%       are needed for the latter.
if second

    optab2 = cell(nt,ntraj);

    % Modifications regarding profiles
    % (which should not include changes to conventional gen, i.e.,
    % gens that are neither wind units nor loads because it could
    % generate conflicts when changing the same unit with more than 1
    % contingency-like changes subsequently in the mpsopfl2 fcn)
    for p = 1:nprof2
        optab2 = profile2contabs(optab2, trajdata(p));
    end


    for t = 1:nt
        for j = 1:ntraj
            md.second.tstep(t).OpCondSched(j).tab = optab2{t,j};
        end
    end

end

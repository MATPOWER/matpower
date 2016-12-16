function [PR_REP, PR_REL, PR_ADD, PR_TCONT, PR_TYPES, PR_TMPCD,...
    PR_TXGD, PR_TCTD, PR_TSTGD, PR_CHGTYPES] = idx_profile
%IDX_PROFILE Defines constants used by Profiles.
%
%	[PR_REP, PR_REL, PR_ADD, PR_TCONT, PR_TYPES, PR_TMPCD,...
%   PR_TXGD, PR_TCTD, PR_TSTGD, PR_CHGTYPES] = idx_profile;
%
%   Indicates and defines numeric and string id's for the types of profiles
%   that can be created to modify input data for MOST across time and
%   scenarios. Some types may require that a table is further specified,
%   i.e., a specific table of data whose content is to be modified by the
%   profile. Rows (rows) and columns (col) fields of each profile further
%   specify which entries of the data in question is to be modified by the
%   profile. The change type field (chgtype) indicates how the change is
%   applied, and finally, the values field (values) contains the new values
%   to be used.
%
%   type:   string field of each profile struct where the id of the type of
%           the profile is stored. Possible values are:
%               string 'mpcData' indicates changes on fields of mpc
%                   (e.g. bus or gen tables)
%               string 'xGenData' indicates changes on xgd fields
%                   (e.g. offers or commitment vars)
%               string 'ContingencyData' indicates changes on the set of
%                   contingencies to be applied 
%               string 'StorageData' indicates changes on the storage struct
% 
%   table:  scalar or string field of each profile struct where the id/name
%           of the table/field that needs to be modified is stored.
%           Possible values are:
%           - if type is 'mpcData' you may modify the following tables
%               (specified with a scalar):
%               scalar CT_TBUS change in bus table
%               scalar CT_TGEN change in gen table
%               scalar CT_TBRCH change in branch table
%               scalar CT_TAREABUS area-wide change in bus table
%               scalar CT_TAREAGEN area-wide change in gen table
%               scalar CT_TAREABRCH area-wide change in branch table
%               scalar CT_TLOAD load modification (bus and/or gen tables)
%               scalar CT_TAREALOAD area-wide change in load (bus/gen tables)
%               scalar CT_TGENCOST change in gencost table
%               scalar CT_TAREAGENCOST area-wide change in gencost table
%           - if type is 'xGenData' you may modify the following fields:
%               string 'CommitSched'
%               string 'InitialPg'
%               string 'RampWearCostCoeff'
%               string 'PositiveActiveReservePrice'
%               string 'PositiveActiveReserveQuantity'
%               string 'NegativeActiveReservePrice'
%               string 'NegativeActiveReserveQuantity'
%               string 'PositiveActiveDeltaPrice'
%               string 'NegativeActiveDeltaPrice'
%               string 'PositiveLoadFollowReservePrice'
%               string 'PositiveLoadFollowReserveQuantity'
%               string 'NegativeLoadFollowReservePrice'
%               string 'NegativeLoadFollowReserveQuantity'
%           - if type is 'ContingencyData' you may modify the following
%             tables:
%               scalar PR_TCONT indicates a modification of the subset of
%               contingencies (of the master table of contingencies) to be
%               applied at each time period and/or scenario.
%           - if type is 'StorageData' you may modify the following fields:
%               string 'MinStorageLevel' 
%               string 'MaxStorageLevel'
%               string 'OutEff'
%               string 'InEff'
%               string 'LossFactor'
%               string 'rho'
%
%   rows:   numeric vector field of each profile struct where the row
%           numbers (1st-dim) of the array to be modified are stored. Row
%           numbers usually represent specifically which gens, branches,
%           buses, contingencies (by the labels), or area will be affected
%           by the modification across time or across scenarios. A value of
%           "0" in this field has the special meaning of "apply to all
%           rows". For an area-wide type of change, the area number is
%           stored here instead.
%
%   col:    scalar field of each profile struct where the id of the
%           parameter to be modified is stored. This id may indicate a
%           column number or some other parameter depending
%           on the type and the table/field. Possible values are:
%           - if type is 'mpcData' the parameters you may modify depend on
%             the table chosen, see IDX_CT on CT_COL for details.
%           - if type is 'xGenData', 'ContingencyData', or 'StorageData',
%             then col is completely ignored.
%
%   chgtype: scalar field of each profile struct where the id of the type of
%           change to be made is stored. Possible values are:
%               scalar PR_REP replaces old values by new 'values'
%               scalar PR_REL multiplies old value by factors in 'vales'
%               scalar PR_ADD adds entries in 'values' field to old value
% 
%   values: numeric array field of each profile struct where the new
%           values/factors, i.e., the profile itself, is stored. This array
%           must have 3 dimensions in a pre-defined order: [nt nj_max n]
%               (i) dimension corresponding to time periods
%              (ii) dimension corresponding to scenarios
%             (iii) dimension corresponding to elements indicated by 'rows'
%           A singleton dimension in 'values' not matching with nt==1,
%           nj_max==1 or length(profile.rows)==1 is interpreted as "apply
%           to all" whenever the parameter being modified allows such an
%           expansion. These "expansions" occur within loadmd().

%   MOST
%   Copyright (c) 2013-2016, Power Systems Engineering Research Center (PSERC)
%   by Daniel Munoz-Alvarez and Ray Zimmerman, PSERC Cornell
%
%   This file is part of MOST.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/most for more info.

%%
[CT_LABEL, CT_PROB, CT_TABLE, CT_TBUS, CT_TGEN, CT_TBRCH, ...
    CT_TAREABUS, CT_TAREAGEN, CT_TAREABRCH, CT_ROW, CT_COL, CT_CHGTYPE, ...
    CT_REP, CT_REL, CT_ADD, CT_NEWVAL, CT_TLOAD, CT_TAREALOAD, ...
    CT_LOAD_ALL_PQ, CT_LOAD_FIX_PQ, CT_LOAD_DIS_PQ, CT_LOAD_ALL_P, ...
    CT_LOAD_FIX_P, CT_LOAD_DIS_P, CT_TGENCOST, CT_TAREAGENCOST, ...
    CT_MODCOST_F, CT_MODCOST_X] = idx_ct;

%% types of profiles

PR_TYPES        = {...                      %% cell array of types
                    'mpcData',...           %% changes to mpc
                    'xGenData',...          %% changes to xgd
                    'ContingencyData',...   %% changes to ct_subset
                    'StorageData',...       %% changes to Storage
                    }; 

%% labels for modifiable tables
% - tables for type 'mpcData'
PR_TMPCD        = [ ... %% vector of tables for 'mpcData'
                    CT_TBUS, ...
                    CT_TGEN, ...
                    CT_TBRCH, ...
                    CT_TAREABUS, ...
                    CT_TAREAGEN, ...
                    CT_TAREABRCH, ...
                    CT_TLOAD, ...
                    CT_TAREALOAD, ...
                    CT_TGENCOST, ...
                    CT_TAREAGENCOST ...
                    ];

% - tables for type xGenData
PR_TXGD         = { ... %% cell array of fields for xGenData
                    'CommitSched',...
                    'InitialPg',...
                    'RampWearCostCoeff',...
                    'PositiveActiveReservePrice',...
                    'PositiveActiveReserveQuantity',...
                    'NegativeActiveReservePrice',...
                    'NegativeActiveReserveQuantity',...
                    'PositiveActiveDeltaPrice',...
                    'NegativeActiveDeltaPrice',...
                    'PositiveLoadFollowReservePrice',...
                    'PositiveLoadFollowReserveQuantity',...
                    'NegativeLoadFollowReservePrice',...
                    'NegativeLoadFollowReserveQuantity',...
                    'CommitKey',...
                    'InitialState',...
                    'MinUp',...
                    'MinDown',...
                    };

% - tables for type 'ContingencyData'
PR_TCONT        = 1;    %% ct_subset of contingencies
PR_TCTD         = [ ... %% vector of tables for 'ContingencyData'
                    PR_TCONT
                    ];

% - tables for type storage
PR_TSTGD        = { ... %% cell array of fields for 'StorageData'
                    'MinStorageLevel',...
                    'MaxStorageLevel',...
                    'OutEff',...
                    'InEff',...
                    'LossFactor',...
                    'rho',...
                    };

%% named values for chgtype entry
PR_REP      = 1;        %% replace old values with new ones in field 'values'
PR_REL      = 2;        %% multiply old values by factors in field 'values'
PR_ADD      = 3;        %% add value in field 'values' to old values
PR_CHGTYPES = [ ...     %% vector of chgtypes
                PR_REP,...
                PR_REL,...
                PR_ADD,...
                ];
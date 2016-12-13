function [CT_LABEL, CT_PROB, CT_TABLE, CT_TBUS, CT_TGEN, CT_TBRCH, ...
    CT_TAREABUS, CT_TAREAGEN, CT_TAREABRCH, CT_ROW, CT_COL, CT_CHGTYPE, ...
    CT_REP, CT_REL, CT_ADD, CT_NEWVAL, CT_TLOAD, CT_TAREALOAD, ...
    CT_LOAD_ALL_PQ, CT_LOAD_FIX_PQ, CT_LOAD_DIS_PQ, CT_LOAD_ALL_P, ...
    CT_LOAD_FIX_P, CT_LOAD_DIS_P, CT_TGENCOST, CT_TAREAGENCOST, ...
    CT_MODCOST_F, CT_MODCOST_X] = idx_ct
%IDX_CT  Defines constants for named column indices to changes table
%
%   [CT_LABEL, CT_PROB, CT_TABLE, CT_TBUS, CT_TGEN, CT_TBRCH, CT_TAREABUS, ...
%   CT_TAREAGEN, CT_TAREABRCH, CT_ROW, CT_COL, CT_CHGTYPE, CT_REP, ...
%   CT_REL, CT_ADD, CT_NEWVAL, CT_TLOAD, CT_TAREALOAD, CT_LOAD_ALL_PQ, ...
%   CT_LOAD_FIX_PQ, CT_LOAD_DIS_PQ, CT_LOAD_ALL_P, CT_LOAD_FIX_P, ...
%   CT_LOAD_DIS_P, CT_TGENCOST, CT_TAREAGENCOST, CT_MODCOST_F, ...
%   CT_MODCOST_X] = idx_ct;
%
% CT_LABEL: column of changes table where the change set label is stored
%
% CT_PROB:  column of changes table where the probability of the
%           change set is stored
%
% CT_TABLE: column of the changes table where the type of system data 
%           table to be modified is stored;
%           type CT_TBUS indicates bus table
%           type CT_TGEN indicates gen table
%           type CT_TBRCH indicates branch table
%           type CT_TLOAD indicates a load modification (bus and/or gen tables)
%           type CT_TAREABUS indicates area-wide change in bus table
%           type CT_TAREAGEN indicates area-wide change in generator table
%           type CT_TAREABRCH indicates area-wide change in branch table
%           type CT_TAREALOAD indicates area-wide change in load
%                             (bus and/or gen tables)
%
% CT_ROW:   column of changes table where the row number in the data
%           table to be modified is stored. A value of "0" in this column
%           has the special meaning "apply to all rows".  For an area-wide
%           type of change, the area number is stored here instead.
%
% CT_COL:   column of changes table where the number of the column in
%           the data table to be modified is stored
%           For CT_TLOAD and CT_TAREALOAD, the value entered in this column
%           is one of the following codes (or its negative), rather than
%           a column index:
%           type CT_LOAD_ALL_PQ modify all loads, real & reactive
%           type CT_LOAD_FIX_PQ modify only fixed loads, real & reactive
%           type CT_LOAD_DIS_PQ modify only dispatchable loads, real & reactive
%           type CT_LOAD_ALL_P modify all loads, real only
%           type CT_LOAD_FIX_P modify only fixed loads, real only
%           type CT_LOAD_DIS_P modify only dispatchable loads, real only
%           If the negative of one of these codes is used, then any affected
%           dispatchable loads will have their costs scaled as well.
%           For CT_TGENCOST and CT_TAREAGENCOST, in addition to an actual
%           column index, this value can also take one of the following
%           codes to indicate a scaling (CT_REL change type) or shifting
%           (CT_ADD change type) of the specified cost functions:
%           type CT_MODCOST_F scales or shifts the cost function vertically
%           type CT_MODCOST_X scales or shifts the cost function horizontally
%           See also MODCOST.
%
% CT_CHGTYPE: column of changes table where the type of change to
%           be made is stored:
%           type CT_REP replaces old value by value in CT_NEWVAL column
%           type CT_REL multiplies old value by factor in CT_NEWVAL column
%           type CT_ADD adds value in CT_NEWVAL column to old value

%   MATPOWER
%   Copyright (c) 2000-2016, Power Systems Engineering Research Center (PSERC)
%   by Carlos E. Murillo-Sanchez, PSERC Cornell & Universidad Nacional de Colombia
%   and Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

%% column labels for changes table
CT_LABEL    = 1;    %% change set label
CT_PROB     = 2;    %% change set probability
CT_TABLE    = 3;    %% type of table to be modified (see possible values below)
CT_ROW      = 4;    %% number of the row to be modified (0 means all rows)
CT_COL      = 5;    %% number of the column to be modified
                    %% (for some values in CT_TABLE column, this can be a
                    %% special code instead of an actual column index)
CT_CHGTYPE  = 6;    %% type of parameter modification to be made
                    %% (see possible values below)
CT_NEWVAL   = 7;    %% quantity to use for replacement value, scale factor
                    %% or shift amount

%% named values for CT_TABLE entry
CT_TBUS         = 1;    %% bus table
CT_TGEN         = 2;    %% gen table
CT_TBRCH        = 3;    %% branch table
CT_TAREABUS     = 4;    %% area-wide change in bus table
CT_TAREAGEN     = 5;    %% area-wide change in gen table
CT_TAREABRCH    = 6;    %% area-wide change in branch table
CT_TLOAD        = 7;    %% single bus load change
CT_TAREALOAD    = 8;    %% area-wide bus load change
CT_TGENCOST     = 9;    %% gencost table
CT_TAREAGENCOST = 10;   %% area-wide change in gencost table

%% named values for CT_CHGTYPE entry
CT_REP      = 1;        %% replace old value with new one in column CT_NEWVAL
CT_REL      = 2;        %% multiply old value by factor in column CT_NEWVAL
CT_ADD      = 3;        %% add value in column CT_NEWVAL to old value

%% codes for CT_COL entry when CT_TABLE entry is CT_TLOAD or CT_TAREALOAD
CT_LOAD_ALL_PQ = 1;     %% all loads, real and reactive
CT_LOAD_FIX_PQ = 2;     %% only fixed loads, real and reactive
CT_LOAD_DIS_PQ = 3;     %% only dispatchable loads, real and reactive
CT_LOAD_ALL_P  = 4;     %% all loads, real only
CT_LOAD_FIX_P  = 5;     %% only fixed loads, real only
CT_LOAD_DIS_P  = 6;     %% only dispatchable loads, real only

%% codes for CT_COL entry when CT_TABLE entry is CT_TGENCOST or CT_TAREAGENCOST
CT_MODCOST_F = -1;      %% scale or shift cost function vertically
CT_MODCOST_X = -2;      %% scale or shift cost function horizontally

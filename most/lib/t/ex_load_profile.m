function loadprofile = ex_load_profile
%EX_LOAD_PROFILE  Example load profile data file for stochastic unit commitment.

%   MOST
%   Copyright (c) 2015-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MOST.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/most for more info.

%% define constants
[CT_LABEL, CT_PROB, CT_TABLE, CT_TBUS, CT_TGEN, CT_TBRCH, CT_TAREABUS, ...
    CT_TAREAGEN, CT_TAREABRCH, CT_ROW, CT_COL, CT_CHGTYPE, CT_REP, ...
    CT_REL, CT_ADD, CT_NEWVAL, CT_TLOAD, CT_TAREALOAD, CT_LOAD_ALL_PQ, ...
    CT_LOAD_FIX_PQ, CT_LOAD_DIS_PQ, CT_LOAD_ALL_P, CT_LOAD_FIX_P, ...
    CT_LOAD_DIS_P, CT_TGENCOST, CT_TAREAGENCOST, CT_MODCOST_F, ...
    CT_MODCOST_X] = idx_ct;

loadprofile = struct( ...
    'type', 'mpcData', ...
    'table', CT_TLOAD, ...
    'rows', 0, ...
    'col', CT_LOAD_ALL_PQ, ...
    'chgtype', CT_REP, ...
    'values', [] );
loadprofile.values(:, 1, 1) = [
	440;
	480;
	540
	525;
	500;
	450;
	400;
	350;
	300;
	325;
	375;
	425;
]-1e-3;

function profiles = uniformwindprofile(nt, nj, n);
%UNIFORMWINDPROFILE Creates a wind profile with evenly spaced capacity values.
%
%   PROFILES = UNIFORMWINDPROFILE(NT, NJ, N)
%
%   Returns a Profile struct to modify the PMAX of a set N generators.
%   The profile has NT periods and NJ scenarios, where the NJ scale
%   factors are evenly spaced from 0 to 1. All inputs are optional.
%   Defaults are NT = 24, NJ = 5, N = 1.
%
%   See also IDX_PROFILE, GET_PROFILES.

%   MOST
%   Copyright (c) 2013-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MOST.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/most for more info.

%% define constants
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[CT_LABEL, CT_PROB, CT_TABLE, CT_TBUS, CT_TGEN, CT_TBRCH, CT_TAREABUS, ...
    CT_TAREAGEN, CT_TAREABRCH, CT_ROW, CT_COL, CT_CHGTYPE, CT_REP, ...
    CT_REL, CT_ADD, CT_NEWVAL, CT_TLOAD, CT_TAREALOAD, CT_LOAD_ALL_PQ, ...
    CT_LOAD_FIX_PQ, CT_LOAD_DIS_PQ, CT_LOAD_ALL_P, CT_LOAD_FIX_P, ...
    CT_LOAD_DIS_P, CT_TGENCOST, CT_TAREAGENCOST, CT_MODCOST_F, ...
    CT_MODCOST_X] = idx_ct;

%% set up default input args
if nargin < 3
    n = [];
    if nargin < 2
        nj = [];
        if nargin < 1
            nt = [];
        end
    end
end
if isempty(nt)
    nt = 24;
end
if isempty(nj)
    nj = 5;
end
if isempty(n)
    n = 1;
end

%% initialize profile
profiles = struct( ...
    'type', 'mpcData', ...
    'table', CT_TGEN, ...
    'rows', 1, ...
    'col', PMAX, ...
    'chgtype', CT_REL, ...
    'values', [] );

c = (0:1/(nj-1):1);
profiles.values = repmat(c, [nt 1 n]);

%% default looks like ...
% profiles.values(:, :, 1) = [
% 	0	0.25	0.5	0.75	1;
% 	0	0.25	0.5	0.75	1;
% 	0	0.25	0.5	0.75	1;
% 	0	0.25	0.5	0.75	1;
% 	0	0.25	0.5	0.75	1;
% 	0	0.25	0.5	0.75	1;
% 	0	0.25	0.5	0.75	1;
% 	0	0.25	0.5	0.75	1;
% 	0	0.25	0.5	0.75	1;
% 	0	0.25	0.5	0.75	1;
% 	0	0.25	0.5	0.75	1;
% 	0	0.25	0.5	0.75	1;
% 	0	0.25	0.5	0.75	1;
% 	0	0.25	0.5	0.75	1;
% 	0	0.25	0.5	0.75	1;
% 	0	0.25	0.5	0.75	1;
% 	0	0.25	0.5	0.75	1;
% 	0	0.25	0.5	0.75	1;
% 	0	0.25	0.5	0.75	1;
% 	0	0.25	0.5	0.75	1;
% 	0	0.25	0.5	0.75	1;
% 	0	0.25	0.5	0.75	1;
% 	0	0.25	0.5	0.75	1;
% 	0	0.25	0.5	0.75	1;
% ];

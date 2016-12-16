function profiles = getprofiles(profilesi, profiles0, idx)
%GETPROFILES  Loads profiles, optionally updating 'rows' via an index mapping.
%
%   PROFILES = GETPROFILES(PROFILESI);
%   PROFILES = GETPROFILES(PROFILESI, PROFILES0);
%   PROFILES = GETPROFILES(PROFILESI, IDX);
%   PROFILES = GETPROFILES(PROFILESI, PROFILES0, IDX);
%
%   Loads a profile or set of profiles from a struct, MAT-file or M-file,
%   optionally using the N-dimensional index vector IDX to modify any
%   non-zero values in ROWS fields so that the corresponding ROWS
%   field in the returned PROFILES is equal to IDX(ROWS). This makes it
%   easy to use profiles defined for a particular set of generators, for
%   example, that are added to a case as a group as in the example below.
%
%   Inputs:
%       PROFILESI : a PROFILE struct or the name of an M-file
%                   or MAT-file that returns one. See IDX_PROFILE for details.
%       PROFILES0 : (optional) a profiles struct to which the newly
%                   loaded profiles specified in PROFILESI will be appended.
%       IDX : (optional) N dimensional index vector used to map non-zero
%             values in ROWS fields.
%
%   Output:
%       PROFILES : Resulting profile struct
%
%   Example: Load a MATPOWER case, add some wind generators, with a
%            profile that applies only to the rows of the added wind
%            units.
%
%       mpc = loadcase('mycase');
%       xgd = loadxgendata('myxgendata');
%       [iwind, mpc, xgd] = addwind('mywindunits', mpc, xgd);
%       profiles = getprofiles('mywindprofile', iwind);
%       md = loadmd(mpc, 'mytransmat', xgd, [], [], profiles);
%
%   See also APPLY_PROFILE, IDX_PROFILE.

%   MOST
%   Copyright (c) 2012-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MOST.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/MATPOWER/most for more info.

%% process input args
if nargin < 2
    profiles0 = [];
    idx = [];
elseif nargin < 3 
    if ~isstruct(profiles0) && ~isempty(profiles0)
        idx = profiles0;
        profiles0 = [];
    else
        idx = [];
    end
end

%% load the profiles
fields = {'type', 'table', 'rows', 'col', 'chgtype', 'values'};
profiles = loadgenericdata(profilesi, 'struct', fields, 'profiles');

%% use idx to map rows
if ~isempty(idx)
    n = length(idx);
    for p = 1:length(profiles)
        k = find(profiles(p).rows);
        if k > 0
            if max(profiles(p).rows(k)) <= n
                profiles(p).rows(k) = idx(profiles(p).rows(k));
            else
                error('getprofiles: to map ROWS using IDX, ROWS must not contain values exceeding LENGTH(IDX) (%d)', n);
            end
        end
    end
end

%% stack profiles
if ~isempty(profiles0)
    profiles = [profiles0; profiles];
end

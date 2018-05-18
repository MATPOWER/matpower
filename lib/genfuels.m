function gf = genfuels()
%GENFUELS   Return list of standard values for generator fuel types
%
%   GF = GENFUELS()
%
%   Returns a cell array of strings containing the following standard
%   generator fuel types for use in the optional MPC.GENFUEL field of the
%   MATPOWER case struct. This is to be considered an unordered list,
%   where the position of a particular fuel type in the list is not
%   defined and is therefore subject to change.
%
%       biomass     - Biomass
%       coal        - Coal
%       dfo         - Distillate Fuel Oil (Diesel, FO1, FO2, FO4)
%       geothermal  - Geothermal
%       hydro       - Hydro
%       hydrops     - Hydro Pumped Storage
%       jetfuel     - Jet Fuel
%       lng         - Liquefied Natural Gas
%       ng          - Natural Gas
%       nuclear     - Nuclear
%       oil         - Unspecified Oil
%       refuse      - Refuse, Municipal Solid Waste
%       rfo         - Residual Fuel Oil (FO5, FO6)
%       solar       - Solar
%       syncgen     - Synchronous Condensor
%       wasteheat   - Waste Heat
%       wind        - Wind
%       wood        - Wood or Wood Waste
%       other       - Other
%       unknown     - Unknown
%       dl          - Dispatchable Load
%       ess         - Energy Storage System
%
%   Example:
%       if ~ismember(mpc.genfuel{k}, genfuels())
%           error('unknown fuel type');
%       end
%
%   See also GENTYPES, SAVECASE.

%   MATPOWER
%   Copyright (c) 2017, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

gf = { ...
    'biomass';
    'coal';
    'dfo';
    'geothermal';
    'hydro';
    'hydrops';
    'jetfuel';
    'lng';
    'ng';
    'nuclear';
    'oil';
    'refuse';
    'rfo';
    'solar';
    'syncgen';
    'wasteheat';
    'wind';
    'wood';
    'other';
    'unknown';
    'dl';
    'ess';
};

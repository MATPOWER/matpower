function om = build_cost_params(om, force)
%BUILD_COST_PARAMS  Builds and saves the full generalized cost parameters.
%
%   -----  DEPRECATED - This method is no longer needed,       -----
%   -----               incorporated into PARAMS_LEGACY_COST   -----
%
%   OM.BUILD_COST_PARAMS()
%   OM.BUILD_COST_PARAMS('force')
%
%   Builds the full set of cost parameters from the individual named
%   sub-sets added via ADD_COSTS. Skips the building process if it has
%   already been done, unless a second input argument is present.
%
%   These cost parameters can be retrieved by calling GET_COST_PARAMS
%   and the user-defined costs evaluated by calling COMPUTE_COST.
%
%   See also OPT_MODEL, ADD_COSTS, PARAMS_LEGACY_COST, GET_COST_PARAMS,
%   COMPUTE_COST.

%   MATPOWER
%   Copyright (c) 2008-2017, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if nargin > 1
    om.cost.params = [];        %% clear cache
end
cp = om.params_legacy_cost();   %% compute and cache full cost parameters

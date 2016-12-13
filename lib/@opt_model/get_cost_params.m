function cp = get_cost_params(om, name, idx)
%GET_COST_PARAMS  Returns the cost parameter struct for user-defined costs.
%   CP = GET_COST_PARAMS(OM)
%   CP = GET_COST_PARAMS(OM, NAME)
%   CP = GET_COST_PARAMS(OM, NAME, IDX)
%
%   Requires calling BUILD_COST_PARAMS first to build the full set of
%   parameters. Returns the full cost parameter struct for all user-defined
%   costs that incorporates all of the named cost sets added via ADD_COSTS,
%   or, if a name is provided it returns the cost struct corresponding to
%   the named set of cost rows (N still has full number of columns).
%
%   The cost parameters are returned in a struct with the following fields:
%       N      - nw x nx sparse matrix
%       Cw     - nw x 1 vector
%       H      - nw x nw sparse matrix (optional, all zeros by default)
%       dd, mm - nw x 1 vectors (optional, all ones by default)
%       rh, kk - nw x 1 vectors (optional, all zeros by default)
%
%   See also OPT_MODEL, ADD_COSTS, BUILD_COST_PARAMS, COMPUTE_COST.

%   MATPOWER
%   Copyright (c) 2008-2016, Power Systems Engineering Research Center (PSERC)
%   by Ray Zimmerman, PSERC Cornell
%
%   This file is part of MATPOWER.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See http://www.pserc.cornell.edu/matpower/ for more info.

if ~isfield(om.cost.params, 'N')
    error('@opt_model/get_cost_params: must call build_cost_params first');
end

cp = om.cost.params;

if nargin > 1
    if getN(om, 'cost', name)
        if nargin < 3 || isempty(idx)
            if prod(size(om.cost.idx.i1.(name))) == 1
                i1 = om.cost.idx.i1.(name);
                iN = om.cost.idx.iN.(name);
            else
                error('@opt_model/get_cost_params: cost set ''%s'' requires an idx arg', name);
            end
        else
            s1 = substruct('.', name, '()', idx);
            i1 = subsref(om.cost.idx.i1, s1);
            iN = subsref(om.cost.idx.iN, s1);
        end
        cp.N  = cp.N(i1:iN,:);
        cp.Cw = cp.Cw(i1:iN);
        cp.H  = cp.H(i1:iN,i1:iN);
        cp.dd = cp.dd(i1:iN);
        cp.rh = cp.rh(i1:iN);
        cp.kk = cp.kk(i1:iN);
        cp.mm = cp.mm(i1:iN);
    end
end
